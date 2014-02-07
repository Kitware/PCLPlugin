/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkAnnotateOBBs.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkAnnotateOBBs.h"

#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkSmartPointer.h"
#include "vtkAppendPolyData.h"
#include "vtkThresholdPoints.h"
#include "vtkOutlineSource.h"
#include "vtkOBBTree.h"

#include <Eigen/Geometry>

//----------------------------------------------------------------------------
namespace {

  static void ComputeOBB(vtkPolyData* polyData, double cornerData[24])
  {
    Eigen::Vector3d origin, x, y, z, sizes;
    vtkSmartPointer<vtkOBBTree> obb = vtkSmartPointer<vtkOBBTree>::New();
    obb->ComputeOBB(polyData->GetPoints(), origin.data(), x.data(), y.data(), z.data(), sizes.data());

    // make right handed
    Eigen::Vector3d rightHandedZ = x.cross(y);
    if (z.dot(rightHandedZ) < 0)
      {
      origin += z;
      z *= -1;
      }

    Eigen::Vector3d corners[8];
    corners[0] = origin;
    corners[1] = origin+x;
    corners[2] = origin+y;
    corners[3] = origin+x+y;
    corners[4] = z+origin;
    corners[5] = z+origin+x;
    corners[6] = z+origin+y;
    corners[7] = z+origin+x+y;
    for (int i = 0; i < 8; ++i)
      {
      cornerData[i*3+0] = corners[i][0];
      cornerData[i*3+1] = corners[i][1];
      cornerData[i*3+2] = corners[i][2];
      }
  }

}


//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkAnnotateOBBs);

//----------------------------------------------------------------------------
vtkAnnotateOBBs::vtkAnnotateOBBs()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  this->AnnotateLabelZero = true;

  this->SetInputArrayToProcess(
    0,0,0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
    vtkDataSetAttributes::SCALARS);
}

//----------------------------------------------------------------------------
vtkAnnotateOBBs::~vtkAnnotateOBBs()
{
}

//----------------------------------------------------------------------------
int vtkAnnotateOBBs::RequestData(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDataArray* labels = this->GetInputArrayToProcess(0, inputVector);
  if (!labels)
    {
    vtkErrorMacro("Could not get input scalars.");
    return 0;
    }

  vtkSmartPointer<vtkPolyData> inputCopy = vtkSmartPointer<vtkPolyData>::New();
  inputCopy->ShallowCopy(input);


  vtkIdType minLabel = labels->GetRange()[0];
  vtkIdType maxLabel = labels->GetRange()[1];

  if (minLabel == 0 && !this->AnnotateLabelZero)
    {
    minLabel = 1;
    }

  vtkNew<vtkAppendPolyData> appendFilter;

  for (vtkIdType i = minLabel; i <= maxLabel; ++i)
    {
    vtkNew<vtkThresholdPoints> threshold;
    threshold->ThresholdBetween(i, i);
    threshold->SetInput(inputCopy);
    threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, labels->GetName());
    threshold->Update();
    vtkPolyData* labelPoints = threshold->GetOutput();

    double cornerData[24];
    ComputeOBB(labelPoints, cornerData);

    vtkSmartPointer<vtkOutlineSource> outline = vtkSmartPointer<vtkOutlineSource>::New();
    outline->SetBoxTypeToOriented();
    outline->SetCorners(cornerData);
    appendFilter->AddInputConnection(outline->GetOutputPort());
    }

  if (appendFilter->GetNumberOfInputConnections(0))
    {
    appendFilter->Update();
    }

  output->ShallowCopy(appendFilter->GetOutput());
  return 1;
}

//----------------------------------------------------------------------------
void vtkAnnotateOBBs::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
