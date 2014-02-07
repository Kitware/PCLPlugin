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
// .NAME vtkParticlePartitionFilter distribute particle datasets in parallel
// .SECTION Description
// vtkParticlePartitionFilter is a parallel load balancing/partitioning 
// filter for particle datasets. It uses the Zoltan library from the Trilinos 
// package to perform the redistribution.

#ifndef __vtkParticlePartitionFilter_h
#define __vtkParticlePartitionFilter_h

#include "vtkDataObjectAlgorithm.h" // superclass
#include "vtkBoundingBox.h"
#include <vector>
#include "zoltan.h"
#include "vtkSmartPointer.h"

class vtkMultiProcessController;
class vtkPoints;
class vtkIdTypeArray;
class vtkBoundsExtentTranslator;

class VTK_EXPORT vtkParticlePartitionFilter : public vtkDataObjectAlgorithm
{
  public:
    static vtkParticlePartitionFilter *New();
    vtkTypeMacro(vtkParticlePartitionFilter,vtkDataObjectAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    // Description:
    // By default this filter uses the global controller,
    // but this method can be used to set another instead.
    virtual void SetController(vtkMultiProcessController*);

    // Description:
    // Specify the name of a scalar array which will be used to fetch
    // the index of each point. This is necessary only if the particles
    // change position (Id order) on each time step. The Id can be used
    // to identify particles at each step and hence track them properly.
    // If this array is NULL, the global point ids are used.  If an Id
    // array cannot otherwise be found, the point index is used as the ID.
    vtkSetStringMacro(IdChannelArray);
    vtkGetStringMacro(IdChannelArray);

    // Description:
    // The thickness of the region between each partition that is used for 
    // ghost cell exchanges. Any particles within this overlap region of another
    // processor will be duplicated on neighbouring processors (possibly multiple times
    // at corner region overlaps)
    vtkSetMacro(GhostCellOverlap, double);
    vtkGetMacro(GhostCellOverlap, double);
    
//BTX
    // Description:
    // Return the Bounding Box for a partition
    vtkBoundingBox *GetPartitionBoundingBox(int partition);

    // Description:
    // Return the Bounding Box for a partition plus the extended region
    // all around the box where ghost cells might be required/present
    vtkBoundingBox *GetPartitionBoundingBoxWithGhostRegion(int partition);
//ETX

  protected:
     vtkParticlePartitionFilter();
    ~vtkParticlePartitionFilter();

    int  GatherDataTypeInfo(vtkDataSet *input);
    bool GatherDataArrayInfo(vtkDataArray *data, int &datatype, std::string &dataname, int &numComponents);

    // Override to specify support for vtkPointSet input type.
    virtual int FillInputPortInformation(int port, vtkInformation* info);

	  // Override to specify different type of output
	  virtual int FillOutputPortInformation(int vtkNotUsed(port), 
		  vtkInformation* info);

    // Description:
    virtual int RequestInformation(vtkInformation*,
                            vtkInformationVector**,
                            vtkInformationVector*);
    // Description:
    virtual int RequestUpdateExtent(vtkInformation*, 
                                   vtkInformationVector**, 
                                   vtkInformationVector*);
    
    // Description:
    virtual int RequestData(vtkInformation*,
                            vtkInformationVector**,
                            vtkInformationVector*);
//BTX
    vtkSmartPointer<vtkIdTypeArray> GenerateGlobalIds(vtkIdType N, const char *idname);

    typedef struct {
      std::vector<ZOLTAN_ID_TYPE> GlobalIds;
      std::vector<ZOLTAN_ID_TYPE> LocalIds;
      std::vector<int> Procs;
    } GhostPartition;

    typedef std::vector< std::vector<vtkIdType> > ListOfVectors;
    void FindOverlappingPoints(vtkPoints *pts, vtkIdTypeArray *IdArray, GhostPartition &ghostinfo);
    vtkBoundingBox             *LocalBox;
    std::vector<vtkBoundingBox> BoxList;
    std::vector<vtkBoundingBox> BoxListWithGhostRegion;
//ETX
    //
    vtkMultiProcessController *Controller;
    //
    int             UpdatePiece;
    int             UpdateNumPieces;
    vtkIdType       NumberOfLocalPoints;
    char           *IdChannelArray;
    double          GhostCellOverlap;
    vtkBoundsExtentTranslator *ExtentTranslator;
    //
  private:
    vtkParticlePartitionFilter(const vtkParticlePartitionFilter&);  // Not implemented.
    void operator=(const vtkParticlePartitionFilter&);  // Not implemented.
};

#endif
