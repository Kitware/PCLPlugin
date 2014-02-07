/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPCLRadiusOutlierRemoval.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPCLRadiusOutlierRemoval.h"
#include "vtkPCLConversions.h"

#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkSmartPointer.h"
#include "vtkNew.h"

#include <pcl/filters/radius_outlier_removal.h>

namespace {

//-----------------------------------------------------------------------------
pcl::IndicesConstPtr ApplyRadiusOutlierRemoval(pcl::PointCloud<pcl::PointXYZ>::ConstPtr cloud,
                                 double searchRadius,
                                 int neighborsInSearchRadius)
{
  if (!cloud || !cloud->points.size())
    {
    return pcl::IndicesConstPtr(new std::vector<int>);
    }

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloudFiltered (new pcl::PointCloud<pcl::PointXYZ>);
  pcl::RadiusOutlierRemoval<pcl::PointXYZ> outrem(true);
  outrem.setInputCloud(cloud);
  outrem.setRadiusSearch(searchRadius);
  outrem.setMinNeighborsInRadius(neighborsInSearchRadius);
  outrem.filter(*cloudFiltered);

  return outrem.getRemovedIndices();
}

}

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPCLRadiusOutlierRemoval);

//----------------------------------------------------------------------------
vtkPCLRadiusOutlierRemoval::vtkPCLRadiusOutlierRemoval()
{
  this->SearchRadius = 0.3;
  this->NeighborsInSearchRadius = 10;
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkPCLRadiusOutlierRemoval::~vtkPCLRadiusOutlierRemoval()
{
}

//----------------------------------------------------------------------------
int vtkPCLRadiusOutlierRemoval::RequestData(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get input and output data objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // perform outlier removal
  pcl::PointIndices::Ptr inlierIndices;
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = vtkPCLConversions::PointCloudFromPolyData(input);
  pcl::IndicesConstPtr outlierIndices = ApplyRadiusOutlierRemoval(cloud,
                                                             this->SearchRadius,
                                                             this->NeighborsInSearchRadius);

  // pass thru input add labels
  vtkSmartPointer<vtkIntArray> labels = vtkPCLConversions::NewLabelsArray(outlierIndices, input->GetNumberOfPoints());
  labels->SetName("is_outlier");
  output->ShallowCopy(input);
  output->GetPointData()->AddArray(labels);

  return 1;
}

//----------------------------------------------------------------------------
void vtkPCLRadiusOutlierRemoval::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
