/*=========================================================================

  Project                 : vtkCSCS
  Module                  : vtkParticlePartitionFilter.h
  Revision of last commit : $Rev: 884 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2010-04-06 12:03:55 +0200 #$

  Copyright (C) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing
  1) This copyright notice appears on all copies of source code
  2) An acknowledgment appears with any substantial usage of the code
  3) If this code is contributed to any other open source project, it
  must not be reformatted such that the indentation, bracketing or
  overall style is modified significantly.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/
#include "vtkToolkits.h"     // For VTK_USE_MPI
//
#ifdef VTK_USE_MPI
  #include "vtkMPI.h"
  #include "vtkMPIController.h"
  #include "vtkMPICommunicator.h"
#endif
#include "vtkMultiProcessController.h"
#include "vtkParticlePartitionFilter.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPolyData.h"
#include "vtkDataSetAttributes.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkTimerLog.h"
#include "vtkIdTypeArray.h"
#include "vtkBoundingBox.h"
#include "vtkMath.h"
//
#include "vtkBoundsExtentTranslator.h"
//
#include <sstream>
//
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <numeric>
//----------------------------------------------------------------------------
#define JB_DEBUG__
#if defined JB_DEBUG__
#define OUTPUTTEXT(a) std::cout <<(a); std::cout.flush();

  #undef vtkDebugMacro
  #define vtkDebugMacro(a)  \
  { \
    if (this->UpdatePiece>=0) { \
      vtkOStreamWrapper::EndlType endl; \
      vtkOStreamWrapper::UseEndl(endl); \
      vtkOStrStreamWrapper vtkmsg; \
      vtkmsg << "P(" << this->UpdatePiece << "): " a << "\n"; \
      OUTPUTTEXT(vtkmsg.str()); \
      vtkmsg.rdbuf()->freeze(0); \
    } \
  }

  #undef  vtkErrorMacro
  #define vtkErrorMacro(a) vtkDebugMacro(a)  
#endif
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkParticlePartitionFilter);
vtkCxxSetObjectMacro(vtkParticlePartitionFilter, Controller, vtkMultiProcessController);

static int pack_count = 0;
static int size_count = 0;
static int unpack_count = 0;
//----------------------------------------------------------------------------
// Zoltan callback interface
//
// Structure to hold mesh data 
//----------------------------------------------------------------------------
typedef struct{
  vtkPointSet        *Input;
  vtkPointSet        *Output;
  vtkIdType           InputNumberOfLocalPoints;
  vtkIdType           OutputNumberOfLocalPoints;
  vtkIdType           OutputNumberOfPointsFinal;
  vtkIdType          *InputGlobalIds;
  void               *InputPointData;   // float/double
  void               *OutputPointData;  // float/double
  vtkPoints          *OutputPoints; 
  int                 NumberOfFields;
  std::vector<void*>  InputArrayPointers;
  std::vector<void*>  OutputArrayPointers;
  std::vector<int>    ArrayTypeSizes;
  int                 TotalSizePerId;
  vtkIdType           OutPointCount;
} PartitionVariables;

//----------------------------------------------------------------------------
// Application defined query functions (Implementation)
//----------------------------------------------------------------------------
static int get_number_of_objects(void *data, int *ierr)
{
  PartitionVariables *mesh= (PartitionVariables *)data;
  *ierr = ZOLTAN_OK;
  return mesh->InputNumberOfLocalPoints;
}
//----------------------------------------------------------------------------
static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
{
  PartitionVariables *mesh = (PartitionVariables*)data;
  //
  // Return the IDs of our objects, but no weights.
  // Zoltan will assume equally weighted objects.
  //
  for (int i=0; i<mesh->InputNumberOfLocalPoints; i++){
    globalID[i] = mesh->InputGlobalIds[i];
    localID[i] = i;
  }
  *ierr = ZOLTAN_OK;
}
//----------------------------------------------------------------------------
static int get_num_geometry(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return 3;
}
//----------------------------------------------------------------------------
template<typename T>
void get_geometry_list(
  void *data, int sizeGID, int sizeLID, int num_obj, 
  ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
  int num_dim, double *geom_vec, int *ierr)
{
  PartitionVariables *mesh = (PartitionVariables*)data;
  for (int i=0;  i < num_obj ; i++){
    geom_vec[3*i]   = ((T*)(mesh->InputPointData))[3*i+0];
    geom_vec[3*i+1] = ((T*)(mesh->InputPointData))[3*i+1];
    geom_vec[3*i+2] = ((T*)(mesh->InputPointData))[3*i+2];
  }
  *ierr = ZOLTAN_OK;
  return;
}
//----------------------------------------------------------------------------
//
/*
A ZOLTAN_OBJ_SIZE_FN query function returns the size (in bytes) of the data buffer 
that is needed to pack all of a single object's data.
 
Function Type: 	ZOLTAN_OBJ_SIZE_FN_TYPE
Arguments: 	
    data 	Pointer to user-defined data.
   num_gid_entries 	The number of array entries used to describe a single global ID.  
      This value is the maximum value over all processors of the parameter NUM_GID_ENTRIES.
   num_lid_entries 	The number of array entries used to describe a single local ID.  
      This value is the maximum value over all processors of the parameter NUM_LID_ENTRIES.
   global_id 	Pointer to the global ID of the object.
   local_id 	Pointer to the local ID of the object.
   ierr 	Error code to be set by function.
Returned Value: 	
    int 	The size (in bytes) of the required data buffer.
*/
//----------------------------------------------------------------------------
template<typename T>
int zoltan_obj_size_func(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr)
{
  size_count++;
  PartitionVariables *mesh = (PartitionVariables*)data;
  *ierr = ZOLTAN_OK;
  return mesh->TotalSizePerId + sizeof(T)*3;
}
//----------------------------------------------------------------------------
template<typename T>
void zoltan_pack_obj_func(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int dest, int size, char *buf, int *ierr)
{
  pack_count++;
  PartitionVariables *mesh = (PartitionVariables*)data;
  vtkIdType GID = *global_id;
  vtkIdType LID = *local_id;
  //
  for (int i=0; i<mesh->NumberOfFields; i++) {
    int asize = mesh->ArrayTypeSizes[i];
    char *dataptr = (char*)(mesh->InputArrayPointers[i]) + asize*LID;
    memcpy(buf, dataptr, asize);
    buf += asize;
  }
  memcpy(buf, &((T*)(mesh->InputPointData))[(*local_id)*3], sizeof(T)*3);  
  *ierr = ZOLTAN_OK;
  return;
}
//----------------------------------------------------------------------------
template<typename T>
void zoltan_unpack_obj_func(void *data, int num_gid_entries,
  ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr)
{
  unpack_count++;
  if (num_gid_entries != 1) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  PartitionVariables *mesh = (PartitionVariables*)data;
  vtkIdType GID = *global_id;
  //
  vtkPointData *inPD  = mesh->Input->GetPointData();
  vtkPointData *outPD = mesh->Output->GetPointData();
  //
  for (int i=0; i<mesh->NumberOfFields; i++) {
    int asize = mesh->ArrayTypeSizes[i];
    char *dataptr = (char*)(mesh->OutputArrayPointers[i]) + asize*(mesh->OutPointCount);
    memcpy(dataptr, buf, asize);
    buf += asize;
  }
  memcpy(&((T*)(mesh->OutputPointData))[mesh->OutPointCount*3], buf, sizeof(T)*3);  
  mesh->OutPointCount++;
  *ierr = ZOLTAN_OK;
  return;
}
//----------------------------------------------------------------------------
/*
data 	            Pointer to user-defined data.
num_gid_entries 	The number of array entries used to describe a single global ID.  This value is the maximum value over all processors of the parameter NUM_GID_ENTRIES.
num_lid_entries 	The number of array entries used to describe a single local ID.  This value is the maximum value over all processors of the parameter NUM_LID_ENTRIES.
num_import 	      The number of objects that will be received by this processor.
import_global_ids An array of num_import global IDs of objects to be received by this processor. This array may be NULL, as the processor does not necessarily need to know which objects it will receive.
import_local_ids 	An array of num_import local IDs of objects to be received by this processor. This array may be NULL, as the processor does not necessarily need to know which objects it will receive.
import_procs 	    An array of size num_import listing the processor IDs of the source processors. This array may be NULL, as the processor does not necessarily need to know which objects is will receive.
import_to_part 	  An array of size num_import listing the parts to which objects will be imported. This array may be NULL, as the processor does not necessarily need to know from which objects it will receive.
num_export 	      The number of objects that will be sent from this processor to other processors.
export_global_ids An array of num_export global IDs of objects to be sent from this processor.
export_local_ids 	An array of num_export local IDs of objects to be sent from this processor.
export_procs 	    An array of size num_export listing the processor IDs of the destination processors.
export_to_part 	  An array of size num_export listing the parts to which objects will be sent.
ierr 	            Error code to be set by function.
*/
template<typename T>
void zoltan_pre_migrate_pp_func(void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  PartitionVariables *mesh = (PartitionVariables*)data;
  // newTotal = original points - sent away + received
  mesh->OutputNumberOfLocalPoints = mesh->InputNumberOfLocalPoints + num_import - num_export;
  mesh->OutputPoints->SetNumberOfPoints(mesh->OutputNumberOfLocalPoints);
  mesh->OutputPointData = mesh->OutputPoints->GetData()->GetVoidPointer(0);
  vtkPointData *inPD  = mesh->Input->GetPointData();
  vtkPointData *outPD = mesh->Output->GetPointData();
  outPD->CopyAllocate(inPD, mesh->OutputNumberOfLocalPoints);
  //
  for (int i=0; i<mesh->NumberOfFields; i++) {
    vtkDataArray *oarray = mesh->Output->GetPointData()->GetArray(i);
    oarray->SetNumberOfTuples(mesh->OutputNumberOfLocalPoints);
    mesh->OutputArrayPointers.push_back(oarray->GetVoidPointer(0));
  }
  std::vector<bool> alive(mesh->InputNumberOfLocalPoints, true);
  for (vtkIdType i=0; i<num_export; i++) {
    alive[export_local_ids[i]] = false;    
  }
  mesh->OutPointCount = 0;
  for (vtkIdType i=0; i<mesh->InputNumberOfLocalPoints; i++) {
    if (alive[i]) {
      outPD->CopyData(inPD, i, mesh->OutPointCount);
      memcpy(&((T*)(mesh->OutputPointData))[mesh->OutPointCount*3], &((T*)(mesh->InputPointData))[i*3], sizeof(T)*3);
      mesh->OutPointCount++;
    }
  }
}
//----------------------------------------------------------------------------
template <typename T>
void zoltan_pre_ghost_migrate_pp_func(void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  PartitionVariables *mesh = (PartitionVariables*)data;
  // resize points to accept ghost cell additions
  mesh->OutputNumberOfPointsFinal = mesh->OutputNumberOfLocalPoints + num_import;
  mesh->OutputPoints->GetData()->Resize(mesh->OutputNumberOfPointsFinal);
  mesh->OutputPoints->SetNumberOfPoints(mesh->OutputNumberOfPointsFinal);
  mesh->OutputPointData = (T*)(mesh->OutputPoints->GetData()->GetVoidPointer(0));
  // copies are now being made from existing data, so use output as new input points
  mesh->InputPointData = mesh->OutputPointData;
  vtkPointData *inPD  = mesh->Output->GetPointData();
  vtkPointData *outPD = mesh->Output->GetPointData();
  //
  // we must resize all the scalar/vector fields for the point data
  // and we set the input pointer to the output ones so that copies use the correct
  // scalar values from their new positions
  mesh->InputArrayPointers.clear();
  mesh->OutputArrayPointers.clear();
  for (int i=0; i<mesh->NumberOfFields; i++) {
    vtkDataArray *oarray = mesh->Output->GetPointData()->GetArray(i);
    oarray->Resize(mesh->OutputNumberOfPointsFinal);
    oarray->SetNumberOfTuples(mesh->OutputNumberOfPointsFinal);
    mesh->InputArrayPointers.push_back(oarray->GetVoidPointer(0));
    mesh->OutputArrayPointers.push_back(oarray->GetVoidPointer(0));
  }
}
//----------------------------------------------------------------------------
// vtkParticlePartitionFilter :: implementation 
//----------------------------------------------------------------------------
vtkParticlePartitionFilter::vtkParticlePartitionFilter()
{
  this->UpdatePiece         = 0;
  this->UpdateNumPieces     = 1;
  this->NumberOfLocalPoints = 0;
  this->Controller          = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  this->IdChannelArray      = NULL;
  this->GhostCellOverlap    = 0.0;
  this->ExtentTranslator    = vtkBoundsExtentTranslator::New();
}

//----------------------------------------------------------------------------
vtkParticlePartitionFilter::~vtkParticlePartitionFilter()
{
  this->SetController(0);
  if (this->IdChannelArray)
    {
    delete [] this->IdChannelArray;
    this->IdChannelArray = NULL;
    }
}

//----------------------------------------------------------------------------
void vtkParticlePartitionFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkParticlePartitionFilter::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // supports any vtkPointSet type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}
//----------------------------------------------------------------------------
int vtkParticlePartitionFilter::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}
//----------------------------------------------------------------------------
vtkBoundingBox *vtkParticlePartitionFilter::GetPartitionBoundingBox(int partition)
{
  if (partition<this->BoxList.size()) {
    return &this->BoxList[partition];
  }
  vtkErrorMacro(<<"Partition not found in Bounding Box list");
  return NULL;
}
//----------------------------------------------------------------------------
vtkBoundingBox *vtkParticlePartitionFilter::GetPartitionBoundingBoxWithGhostRegion(int partition)
{
  if (partition<this->BoxListWithGhostRegion.size()) {
    return &this->BoxListWithGhostRegion[partition];
  }
  vtkErrorMacro(<<"Partition not found in Bounding Box list");
  return NULL;
}
//----------------------------------------------------------------------------
void vtkParticlePartitionFilter::FindOverlappingPoints(vtkPoints *pts, vtkIdTypeArray *IdArray, GhostPartition &ghostinfo)
{
  vtkIdType N = pts->GetNumberOfPoints();
  for (vtkIdType i=0; i<N; i++) {
    double *pt = pts->GetPoint(i);
    int proc = 0;
    for (std::vector<vtkBoundingBox>::iterator it=this->BoxListWithGhostRegion.begin(); 
      it!=this->BoxListWithGhostRegion.end(); ++it, ++proc) 
    {
      vtkBoundingBox &b = *it;
      if (&b!=this->LocalBox && b.ContainsPoint(pt)) {
        ghostinfo.GlobalIds.push_back(IdArray->GetValue(i));
        ghostinfo.LocalIds.push_back(i);
        ghostinfo.Procs.push_back(proc);
      }
    }
  }
}
//----------------------------------------------------------------------------
vtkSmartPointer<vtkIdTypeArray> vtkParticlePartitionFilter::GenerateGlobalIds(vtkIdType N, const char *idname)
{
  vtkSmartPointer<vtkIdTypeArray> Ids = vtkSmartPointer<vtkIdTypeArray>::New();

  vtkstd::vector<int>       PartialSum(this->UpdateNumPieces+1);
  vtkstd::vector<vtkIdType> PointsPerProcess(this->UpdateNumPieces);
  //
  vtkMPICommunicator* com = vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
  com->AllGather(&N, &PointsPerProcess[0], 1);
  vtkstd::partial_sum(PointsPerProcess.begin(), PointsPerProcess.end(), PartialSum.begin()+1);

  vtkIdType initialValue = PartialSum[this->UpdatePiece];
//  vtkDebugMacro(<< "Id filter rank " << this->UpdatePiece << " Using offset " << initialValue);
//  for (int i=0; i<PartialSum.size(); i++) { vtkDebugMacro(<<PartialSum[i] << " " ; } vtkDebugMacro(<<std::endl;
  //
  Ids->SetNumberOfValues(N);
  for (vtkIdType id=0; id<N; id++) {
    Ids->SetValue(id, id+initialValue);
  }
  Ids->SetName(idname);
  return Ids;
}
//----------------------------------------------------------------------------
struct vtkPPF_datainfo {
  int  datatype;
  int  numC;
  char name[64];
  vtkPPF_datainfo() : datatype(-1), numC(-1) {};
};
//----------------------------------------------------------------------------
bool vtkParticlePartitionFilter::GatherDataArrayInfo(vtkDataArray *data, 
  int &datatype, std::string &dataname, int &numComponents)
{
#ifdef VTK_USE_MPI
  std::vector< vtkPPF_datainfo > datatypes(this->UpdateNumPieces);
  if (data) {
    ((vtkPPF_datainfo*)&datatypes[this->UpdatePiece])->datatype = data->GetDataType();
    ((vtkPPF_datainfo*)&datatypes[this->UpdatePiece])->numC     = data->GetNumberOfComponents();
    strncpy(((vtkPPF_datainfo*)&datatypes[this->UpdatePiece])->name, data->GetName(), 64);
  }
  vtkMPICommunicator* com = vtkMPICommunicator::SafeDownCast(
    this->Controller->GetCommunicator()); 
  int result = com->AllGather((char*)MPI_IN_PLACE, (char*)&datatypes[0], sizeof(vtkPPF_datainfo));
  for (int i=0; i<this->UpdateNumPieces; i++) {
    vtkPPF_datainfo &newdata = datatypes[i];
    if (newdata.datatype!=-1) {
      datatype = newdata.datatype;
      numComponents = newdata.numC;
      dataname = newdata.name;
    }
  }
  return (result == 1) ;
#else
  return 1;
#endif
}
//----------------------------------------------------------------------------
int vtkParticlePartitionFilter::GatherDataTypeInfo(vtkDataSet *input)
{
  vtkPointSet *pInput = vtkPointSet::SafeDownCast(input);
#ifdef VTK_USE_MPI
  if (this->UpdateNumPieces==1) {
      return pInput->GetPoints()->GetDataType();
  }
  std::vector< int > datatypes(this->UpdateNumPieces, -1);
  int datatype = -1;
  if (pInput->GetPoints()) {
    datatypes[this->UpdatePiece] = pInput->GetPoints()->GetDataType();
  }
  vtkMPICommunicator* com = vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator()); 
  int result = com->AllGather((int*)MPI_IN_PLACE, (int*)&datatypes[0], 1);
  for (int i=0; i<this->UpdateNumPieces; i++) {
    int &newdatatype = datatypes[i];
    if (datatype==-1 && newdatatype!=-1) {
      datatype = newdatatype;
    }
    else if (datatype!=-1 && newdatatype!=-1 && newdatatype!=datatype) {
      vtkErrorMacro(<<"Fatal datatype error in Point DataType Gather");
    }
  }
  return datatype;
#else
  return pInput->GetPoints()->GetDataType();
#endif
}
//-------------------------------------------------------------------------
int vtkParticlePartitionFilter::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  int piece, numPieces;
  piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  numPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), piece);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), numPieces);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);

  return 1;
}
//----------------------------------------------------------------------------
int vtkParticlePartitionFilter::RequestInformation(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
#ifdef VTK_USE_MPI
  vtkMPICommunicator *communicator = 
    vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
  int maxpieces = communicator ? communicator->GetNumberOfProcesses() : 1;
#else
  int maxpieces = 1;
#endif

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR(), this->ExtentTranslator);
  //
//  outInfo->Set(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR(),
//               inInfo->Get(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR()));
  //
  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
  //
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()),
               6);
  return 1;
}
//----------------------------------------------------------------------------
int vtkParticlePartitionFilter::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  pack_count = 0;
  size_count = 0;
  unpack_count = 0;
  //
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPointSet     *output = vtkPointSet::GetData(outputVector,0);
  vtkInformation  *inInfo = inputVector[0]->GetInformationObject(0);
  vtkPointSet      *input = vtkPointSet::GetData(inputVector[0]);
  //
  this->UpdatePiece     = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  int ghostLevel        = inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
  vtkDebugMacro(<<"Partition filter " << this->UpdatePiece << " Ghost level " << ghostLevel);
  //
#ifdef VTK_USE_MPI
  vtkMPICommunicator *communicator = vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
  MPI_Comm mpiComm = MPI_COMM_NULL;
  if (communicator) {
    mpiComm = *(communicator->GetMPIComm()->GetHandle());
  }
#else
  int mpiComm = 0;
#endif

  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();

  // Get input
  vtkIdType       numPoints = input->GetNumberOfPoints();
  vtkDataArray    *inPoints = numPoints>0 ? input->GetPoints()->GetData() : NULL;
  vtkDebugMacro(<<"Partitioning on " << this->UpdatePiece << " Points Input : " << numPoints);

  // Setup output
  vtkSmartPointer<vtkPoints>   outPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
  // if input had 0 points, make sure output is still setup correctly (float/double?)
  int pointstype = this->GatherDataTypeInfo(input);
  outPoints->SetDataType(pointstype);
  output->SetPoints(outPoints);

  //
  // We'd like to clamp bounding boxes of all the generated partitions to the original data
  // rather than DBL_MAX/MIN etc used by zoltan (infinite extents)
  //
  double bounds[6];
  input->GetBounds(bounds);
  this->BoxList.clear();
  this->BoxListWithGhostRegion.clear();
  if (this->UpdateNumPieces==1) {
    output->ShallowCopy(input);
    vtkBoundingBox box(bounds);
    this->BoxList.push_back(box);
    // we add a ghost cell region to our boxes
    box.Inflate(this->GhostCellOverlap);
    this->BoxListWithGhostRegion.push_back(box);
    this->ExtentTranslator->SetNumberOfPieces(1);
    // Copy the bounds to our piece to bounds translator
    this->ExtentTranslator->SetBoundsForPiece(0, bounds);
    this->ExtentTranslator->InitWholeBounds();
    //
    return 1;
  }

  double bmin[3] = {bounds[0], bounds[2], bounds[4]};
  double bmax[3] = {bounds[1], bounds[3], bounds[5]};
  if (!vtkMath::AreBoundsInitialized(bounds)) {
    bmin[0] = bmin[1] = bmin[2] = VTK_DOUBLE_MAX;
    bmax[0] = bmax[1] = bmax[2] = VTK_DOUBLE_MIN;
  }
  if (this->Controller) {
    double globalMins[3], globalMaxes[3];
    this->Controller->AllReduce(bmin, globalMins,  3, vtkCommunicator::MIN_OP);
    this->Controller->AllReduce(bmax, globalMaxes, 3, vtkCommunicator::MAX_OP);
    bmin[0] = globalMins[0];  bmax[0] = globalMaxes[0];
    bmin[1] = globalMins[1];  bmax[1] = globalMaxes[1];
    bmin[2] = globalMins[2];  bmax[2] = globalMaxes[2];
  }
  
  //
  // we make a temp copy of the input so we can add Ids if necessary
  //
  vtkSmartPointer<vtkPointSet> inputCopy = input->NewInstance();
  inputCopy->ShallowCopy(input);

  //--------------------------------------------------------------
  // Use Zoltan library to re-partition the particles in parallel
  //--------------------------------------------------------------
  struct Zoltan_Struct *zz;
  int changes, numGidEntries, numLidEntries, numImport=0, numExport=0;
  ZOLTAN_ID_PTR importGlobalGids = NULL;
  ZOLTAN_ID_PTR importLocalGids  = NULL; 
  ZOLTAN_ID_PTR exportGlobalGids = NULL;
  ZOLTAN_ID_PTR exportLocalGids  = NULL;
  int *importProcs  = NULL;
  int *importToPart = NULL;
  int *exportProcs  = NULL;
  int *exportToPart = NULL;
  PartitionVariables mesh;

  float ver;
  int rc = Zoltan_Initialize(0, NULL, &ver);
  if (rc != ZOLTAN_OK){
    printf("Zoltan initialization failed ...\n");
    return 0;
  }

  mesh.Input                    = inputCopy;
  mesh.Output                   = output;
  mesh.InputNumberOfLocalPoints = numPoints;
  mesh.InputPointData           = inPoints ? inPoints->GetVoidPointer(0) : NULL;
  mesh.OutputPoints             = outPoints;
  mesh.TotalSizePerId           = 0;
  mesh.OutPointCount            = 0;

  //
  // if a process has zero points, we need to make dummy point data arrays to allow 
  // space for when data gets sent in from other processes in the zoltan unpack function 
  //
  int NumberOfFields = inputCopy->GetPointData()->GetNumberOfArrays();
  this->Controller->AllReduce(&NumberOfFields, &mesh.NumberOfFields, 1, vtkCommunicator::MAX_OP);
  for (int i=0; i<mesh.NumberOfFields; i++) {
    vtkSmartPointer<vtkDataArray> darray = inputCopy->GetPointData()->GetArray(i);
    //
    int correctType = -1, numComponents = -1;
    std::string correctName;
    this->GatherDataArrayInfo(darray, correctType, correctName, numComponents);
    if (!darray) {
      vtkDebugMacro(<<"NULL data found, used MPI_Gather to find :" 
        << " DataType " << correctType
        << " Name " << correctName.c_str()
        << " NumComponents " << numComponents);
      darray.TakeReference(vtkDataArray::CreateDataArray(correctType));
      darray->SetNumberOfComponents(numComponents);
      darray->SetName(correctName.c_str());
      inputCopy->GetPointData()->AddArray(darray);
    }
  }

  //
  // Global Ids : always do them after other point arrays 
  //
  std::string IdsName;
  if (this->IdChannelArray) {
    IdsName = this->IdChannelArray;
  }
  if (IdsName.empty() || IdsName==std::string("Not available")) {
    IdsName = "PPF_PointIds";
  } 

  vtkSmartPointer<vtkDataArray> Ids = NULL;
  Ids = input->GetPointData()->GetArray(IdsName.c_str());
  if (!Ids) {
    // Try loading the global ids.
    Ids = input->GetPointData()->GetGlobalIds();
  }
  if (!Ids) {
    // Generate our own since none exist
    Ids = this->GenerateGlobalIds(numPoints, IdsName.c_str());
    inputCopy->GetPointData()->AddArray(Ids);
    // and increment the mesh field count
    mesh.NumberOfFields++;
  }
  mesh.InputGlobalIds = vtkIdTypeArray::SafeDownCast(Ids)->GetPointer(0);

  for (int i=0; i<mesh.NumberOfFields; i++) {
    vtkSmartPointer<vtkDataArray> darray = inputCopy->GetPointData()->GetArray(i);
    mesh.InputArrayPointers.push_back(darray->GetVoidPointer(0));
    mesh.ArrayTypeSizes.push_back(darray->GetDataTypeSize());
    mesh.TotalSizePerId += darray->GetDataTypeSize();
  }

  //***************************************************************
  //* Create a Zoltan library structure for this instance of load
  //* balancing.  Set the parameters and query functions that will
  //* govern the library's calculation.  See the Zoltan User's
  //* Guide for the definition of these and many other parameters.
  //***************************************************************

  zz = Zoltan_Create(mpiComm); 

  // we don't need any debug info
  Zoltan_Set_Param(zz, "RCB_OUTPUT_LEVEL", "0");
  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");

  // Method for subdivision
  Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");
  Zoltan_Set_Param(zz, "LB_METHOD",   "RCB");
  //  Zoltan_Set_Param(zz, "LB_METHOD", "PARMETIS");

  // Global and local Ids are a single integer
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");

  // divide into N global and M local partitions
  std::stringstream global;
  global << this->UpdateNumPieces << ends;
  std::stringstream local;
  local << 1 << ends;

  Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", global.str().c_str());
//  Zoltan_Set_Param(zz, "NUM_LOCAL_PARTS",  local.str().c_str());

  // All points have the same weight
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

  // RCB parameters
  //  Zoltan_Set_Param(zz, "PARMETIS_METHOD", "PARTKWAY");
  Zoltan_Set_Param(zz, "RCB_RECOMPUTE_BOX", "0");
  Zoltan_Set_Param(zz, "REDUCE_DIMENSIONS", "0");
  Zoltan_Set_Param(zz, "RCB_MAX_ASPECT_RATIO", "10");

  // we need the cuts to get BBoxes for partitions later
  Zoltan_Set_Param(zz, "KEEP_CUTS", "1");

  // don't allow points on cut to be in different partitions
  // not likely/useful for particle data anyway
  Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "1"); 

  // Let Zoltan do the load balance step automatically
  // particles will be transferred as required between processes
  Zoltan_Set_Param(zz, "AUTO_MIGRATE", "1");  
  
  //
  // Query functions, to provide geometry to Zoltan 
  //
  Zoltan_Set_Num_Obj_Fn(zz,    get_number_of_objects, &mesh);
  Zoltan_Set_Obj_List_Fn(zz,   get_object_list,       &mesh);
  Zoltan_Set_Num_Geom_Fn(zz,   get_num_geometry,      &mesh);
  if (pointstype==VTK_FLOAT) {
    vtkDebugMacro(<<"Using float data pointers ");
    Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list<float>, &mesh);
  }
  else if (pointstype==VTK_DOUBLE) {
    vtkDebugMacro(<<"Using double data pointers ");
    Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list<double>, &mesh);
  }

  //
  // Register functions for packing and unpacking data
  // by migration tools.  
  // GCC has trouble resolving the templated function pointers, so we explicitly
  // declare the types ant then cast them as args
  //    

  
  typedef int  (*zsize_fn) (void *, int , int , ZOLTAN_ID_PTR , ZOLTAN_ID_PTR , int *);
  typedef void (*zpack_fn) (void *, int , int , ZOLTAN_ID_PTR , ZOLTAN_ID_PTR , int , int , char *, int *);
  typedef void (*zupack_fn)(void *, int , ZOLTAN_ID_PTR , int , char *, int *);
  typedef void (*zprem_fn) (void *, int , int , int , ZOLTAN_ID_PTR , ZOLTAN_ID_PTR , int *, int *, int , ZOLTAN_ID_PTR ,
    ZOLTAN_ID_PTR , int *, int *, int *);

  if (pointstype==VTK_FLOAT) {
    zsize_fn  f1 = zoltan_obj_size_func<float>;
    zpack_fn  f2 = zoltan_pack_obj_func<float>;
    zupack_fn f3 = zoltan_unpack_obj_func<float>;
    zprem_fn  f4 = zoltan_pre_migrate_pp_func<float>; 
    Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_FN_TYPE,       (void (*)()) f1, &mesh); 
    Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_FN_TYPE,       (void (*)()) f2, &mesh); 
    Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_FN_TYPE,     (void (*)()) f3, &mesh); 
    Zoltan_Set_Fn(zz, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &mesh); 
  }
  else if (pointstype==VTK_DOUBLE) {
    zsize_fn  f1 = zoltan_obj_size_func<double>;
    zpack_fn  f2 = zoltan_pack_obj_func<double>;
    zupack_fn f3 = zoltan_unpack_obj_func<double>;
    zprem_fn  f4 = zoltan_pre_migrate_pp_func<double>;
    Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_FN_TYPE,       (void (*)()) f1, &mesh);
    Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_FN_TYPE,       (void (*)()) f2, &mesh);
    Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_FN_TYPE,     (void (*)()) f3, &mesh);
    Zoltan_Set_Fn(zz, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &mesh);
  }
  //
  // Zoltan can now partition our particles. 
  // After this returns, we have redistributed particles and the Output holds
  // the list of correct points/fields etc for each process
  //
  rc = Zoltan_LB_Partition(zz, // input (all remaining fields are output)
        &changes,              // 1 if partitioning was changed, 0 otherwise 
        &numGidEntries,        // Number of integers used for a global ID
        &numLidEntries,        // Number of integers used for a local ID
        &numImport,            // Number of vertices to be sent to me
        &importGlobalGids,     // Global IDs of vertices to be sent to me
        &importLocalGids,      // Local IDs of vertices to be sent to me
        &importProcs,          // Process rank for source of each incoming vertex
        &importToPart,         // New partition for each incoming vertex
        &numExport,            // Number of vertices I must send to other processes*/
        &exportGlobalGids,     // Global IDs of the vertices I must send
        &exportLocalGids,      // Local IDs of the vertices I must send
        &exportProcs,          // Process to which I send each of the vertices
        &exportToPart);        // Partition to which each vertex will belong

  if (rc != ZOLTAN_OK){
    printf("Zoltan_LB_Partition NOT OK...\n");
    MPI_Finalize();
    Zoltan_Destroy(&zz);
    exit(0);
  }

  vtkDebugMacro(<<"Partitioning complete on " << this->UpdatePiece << 
    " pack_count : " << pack_count <<
    " size_count : " << size_count <<
    " unpack_count : " << unpack_count 
   );

  MPI_Barrier(mpiComm);
  Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                      &importProcs, &importToPart);
  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                      &exportProcs, &exportToPart);

  //
  // For ghost cells we would like the bounding boxes of each partition
  //
  this->ExtentTranslator->SetNumberOfPieces(this->UpdateNumPieces);
  for (int p=0; p<this->UpdateNumPieces; p++) {
    double bounds[6];
    int ndim;
    if (ZOLTAN_OK==Zoltan_RCB_Box(zz, p, &ndim, &bounds[0], &bounds[2], &bounds[4], &bounds[1], &bounds[3], &bounds[5])) {
      if (bounds[0]==-DBL_MAX) { bounds[0] = bmin[0]; }
      if (bounds[1]== DBL_MAX) { bounds[1] = bmax[0]; }
      if (bounds[2]==-DBL_MAX) { bounds[2] = bmin[1]; }
      if (bounds[3]== DBL_MAX) { bounds[3] = bmax[1]; }
      if (bounds[4]==-DBL_MAX) { bounds[4] = bmin[2]; }
      if (bounds[5]== DBL_MAX) { bounds[5] = bmax[2]; }
      vtkBoundingBox box(bounds);
      this->BoxList.push_back(box);
      //
      // Copy the bounds to our piece to bounds translator
      this->ExtentTranslator->SetBoundsForPiece(p, bounds);
      //
      // Add a ghost cell region to our boxes
      box.Inflate(this->GhostCellOverlap);
      this->BoxListWithGhostRegion.push_back(box);
    }
  }
  this->ExtentTranslator->InitWholeBounds();
  this->LocalBox = &this->BoxListWithGhostRegion[this->UpdatePiece];

  //
  // Find points which overlap other processes' ghost regions
  // note that we must use the 'new' migrated points which are not the same
  // as the original input points (might be bigger/smaller), so get the new IdArray 
  //
  vtkIdTypeArray *newIds = vtkIdTypeArray::SafeDownCast(
    mesh.Output->GetPointData()->GetArray(IdsName.c_str()));
  if (!newIds || newIds->GetNumberOfTuples()!=mesh.OutputPoints->GetNumberOfPoints()) {
    vtkErrorMacro(<<"Fatal : Ids on migrated data corrupted");
    return 0;
  }

  GhostPartition GhostIds;
  this->FindOverlappingPoints(mesh.OutputPoints, newIds, GhostIds);

  //
  // Pass the lists of ghost cells to zoltan so that it
  // can build a list of lists for exchanges between processes
  //
  size_t        num_known = GhostIds.GlobalIds.size(); 
  int           num_found = 0;
  ZOLTAN_ID_PTR found_global_ids = NULL;
  ZOLTAN_ID_PTR found_local_ids  = NULL;
  int          *found_procs      = NULL;
  int          *found_to_part    = NULL;
  //
  rc = Zoltan_Invert_Lists(zz, 
        (int)num_known,
        num_known>0 ? &GhostIds.GlobalIds[0] : NULL,
        num_known>0 ? &GhostIds.LocalIds[0]  : NULL,
        num_known>0 ? &GhostIds.Procs[0]     : NULL,
        num_known>0 ? &GhostIds.Procs[0]     : NULL,
        &num_found,
        &found_global_ids,
        &found_local_ids,
        &found_procs,
        &found_to_part); 

  if (rc != ZOLTAN_OK){
    printf("Zoltan_LB_Partition NOT OK...\n");
    MPI_Finalize();
    Zoltan_Destroy(&zz);
    exit(0);
  }

  //
  // Before sending, we need to change the pre-migrate function as we are now adding
  // extra ghost cells and not starting our lists from a clean slate.
  //
  if (pointstype==VTK_FLOAT) {
    zprem_fn f4 = zoltan_pre_ghost_migrate_pp_func<float>;
    Zoltan_Set_Fn(zz, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &mesh);
  }
  else if (pointstype==VTK_DOUBLE) {
    zprem_fn f4 = zoltan_pre_ghost_migrate_pp_func<double>;
    Zoltan_Set_Fn(zz, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &mesh);
  }

  //
  // Now we can actually send ghost particles between processes
  //
	rc = Zoltan_Migrate (zz,
        (int)num_found,
        found_global_ids,
        found_local_ids,
        found_procs,
        found_to_part,
        (int)num_known,
        num_known>0 ? &GhostIds.GlobalIds[0] : NULL,
        num_known>0 ? &GhostIds.LocalIds[0]  : NULL,
        num_known>0 ? &GhostIds.Procs[0]     : NULL,
        num_known>0 ? &GhostIds.Procs[0]     : NULL
      );

  //
  // Release the arrays allocated during Zoltan_Invert_Lists
  //
  Zoltan_LB_Free_Part(&found_global_ids, &found_local_ids, 
                      &found_procs, &found_to_part);

  //
  // Ghost information : Paraview doesn't let us visualize an array called vtkGhostLevels
  // because it's an 'internal' array, so we make an extra one for debug purposes
  //
  vtkSmartPointer<vtkUnsignedCharArray> GhostArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
  vtkSmartPointer<vtkIntArray> GhostArray2 = vtkSmartPointer<vtkIntArray>::New();
  GhostArray->SetName("vtkGhostLevels");
  GhostArray->SetNumberOfComponents(1);
  GhostArray->SetNumberOfTuples(mesh.OutputNumberOfLocalPoints + num_found);
  GhostArray2->SetName("GhostLevels");
  GhostArray2->SetNumberOfComponents(1);
  GhostArray2->SetNumberOfTuples(mesh.OutputNumberOfLocalPoints + num_found);
  unsigned char *ghost = GhostArray->GetPointer(0);
  int          *ghost2 = GhostArray2->GetPointer(0);
  for (vtkIdType i=0; i<mesh.OutputNumberOfLocalPoints + num_found; i++) {
    if (i<mesh.OutputNumberOfLocalPoints) {
      ghost[i]  = 0;
      ghost2[i] = 0;
    }
    else {
      ghost[i]  = 1;
      ghost2[i] = 1;
    }
  }
  output->GetPointData()->AddArray(GhostArray);
  output->GetPointData()->AddArray(GhostArray2);
  
  //
  //
  //
  vtkDebugMacro(<<"Process " << this->UpdatePiece << " Points Output : " << mesh.OutPointCount);
//  for (int i=0; i<mesh.NumberOfFields; i++) {
//    vtkDataArray *darray = output->GetPointData()->GetArray(i);
//    vtkDebugMacro(<<"Process " << this->UpdatePiece << " Array Output : " << darray->GetNumberOfTuples());
//  }

  //
  // If polydata create Vertices for each point
  //
  if (vtkPolyData::SafeDownCast(output)) {
    vtkIdType *arraydata = vertices->WritePointer(mesh.OutPointCount, 2*mesh.OutPointCount);
    for (int i=0; i<mesh.OutPointCount; i++) {
      arraydata[i*2]   = 1;
      arraydata[i*2+1] = i;
    }
    vtkPolyData::SafeDownCast(output)->SetVerts(vertices);
  }
  //
  //*****************************************************************
  // Free the arrays allocated by Zoltan_LB_Partition, and free
  // the storage allocated for the Zoltan structure.
  //*****************************************************************
  //
  Zoltan_Destroy(&zz);

  this->Controller->Barrier();
  timer->StopTimer();
  vtkDebugMacro(<<"Particle partitioning : " << timer->GetElapsedTime() << " seconds");
  return 1;
}
