/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPCLNormalEstimation.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPCLNormalEstimation.h"
#include "vtkPCLConversions.h"

#include "vtkPolyData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkSmartPointer.h"
#include "vtkNew.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"

#include <pcl/features/normal_3d.h>

//----------------------------------------------------------------------------
namespace {

pcl::PointCloud<pcl::Normal>::Ptr ComputeNormalEstimation(
      pcl::PointCloud<pcl::PointXYZ>::ConstPtr cloud,
      pcl::PointCloud<pcl::PointXYZ>::ConstPtr searchCloud,
      double searchRadius)
{
  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
  pcl::PointCloud<pcl::Normal>::Ptr cloud_normals(new pcl::PointCloud<pcl::Normal>);
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ> ());

  ne.setSearchMethod(tree);
  ne.setInputCloud(cloud);
  if (searchCloud)
    {
    ne.setSearchSurface(searchCloud);
    }

  ne.setRadiusSearch(searchRadius);
  ne.compute(*cloud_normals);
  return cloud_normals;
}

}

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPCLNormalEstimation);

//----------------------------------------------------------------------------
vtkPCLNormalEstimation::vtkPCLNormalEstimation()
{
  this->SearchRadius = 0.1;
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkPCLNormalEstimation::~vtkPCLNormalEstimation()
{
}

//----------------------------------------------------------------------------
int vtkPCLNormalEstimation::FillInputPortInformation(int port, vtkInformation* info)
{
  if (port == 1)
    {
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }

  return this->Superclass::FillInputPortInformation(port, info);
}

//----------------------------------------------------------------------------
int vtkPCLNormalEstimation::RequestData(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));


  vtkInformation *searchInfo = inputVector[1]->GetInformationObject(0);
  vtkPolyData *searchCloudPolyData = 0;
  if (searchInfo)
    {
    searchCloudPolyData = vtkPolyData::SafeDownCast(searchInfo->Get(vtkDataObject::DATA_OBJECT()));
    }

  // early exit if input data has no points
  if (!input->GetNumberOfPoints())
    {
    return 1;  
    }

  // create new normals array
  vtkSmartPointer<vtkFloatArray> normals = vtkSmartPointer<vtkFloatArray>::New();
  normals->SetNumberOfComponents(3);
  normals->SetNumberOfTuples(input->GetNumberOfPoints());
  normals->SetName("normals");

  // pass input thru to output and add new array
  output->ShallowCopy(input);
  output->GetPointData()->AddArray(normals);


  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = vtkPCLConversions::PointCloudFromPolyData(input);
  pcl::PointCloud<pcl::PointXYZ>::Ptr searchCloud;
  if (searchCloudPolyData)
    {
    searchCloud = vtkPCLConversions::PointCloudFromPolyData(searchCloudPolyData);
    }

  pcl::PointCloud<pcl::Normal>::Ptr cloudNormals = ComputeNormalEstimation(cloud,
    searchCloud, this->SearchRadius);

  assert(cloudNormals);
  assert(cloudNormals->size() == normals->GetNumberOfTuples());

  for (size_t i = 0; i < cloudNormals->size(); ++i)
    {
    normals->SetTuple(i, cloudNormals->points[i].normal);
    }

  return 1;
}

//----------------------------------------------------------------------------
void vtkPCLNormalEstimation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
