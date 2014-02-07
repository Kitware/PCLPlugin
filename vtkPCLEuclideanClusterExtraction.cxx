/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPCLEuclideanClusterExtraction.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPCLEuclideanClusterExtraction.h"
#include "vtkPCLConversions.h"

#include "vtkPolyData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPointData.h"

#include <pcl/kdtree/kdtree.h>
#include <pcl/segmentation/extract_clusters.h>

//----------------------------------------------------------------------------
namespace {

void ApplyEuclideanClusterExtraction(
                      pcl::PointCloud<pcl::PointXYZ>::ConstPtr cloud,
                      double clusterTolerance, double minClusterSize, double maxClusterSize,
                      std::vector<pcl::PointIndices>& clusterIndices)
{
  if (!cloud || !cloud->points.size())
    {
    return;
    }

  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
  tree->setInputCloud(cloud);
  pcl::EuclideanClusterExtraction<pcl::PointXYZ> clusterFilter;
  clusterFilter.setClusterTolerance(clusterTolerance);
  clusterFilter.setMinClusterSize(minClusterSize);
  clusterFilter.setMaxClusterSize(maxClusterSize);
  clusterFilter.setSearchMethod(tree);
  clusterFilter.setInputCloud(cloud);
  clusterFilter.extract(clusterIndices);
}

}

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPCLEuclideanClusterExtraction);

//----------------------------------------------------------------------------
vtkPCLEuclideanClusterExtraction::vtkPCLEuclideanClusterExtraction()
{
  this->ClusterTolerance = 0.05;
  this->MinClusterSize = 100;
  this->MaxClusterSize = 100000;
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkPCLEuclideanClusterExtraction::~vtkPCLEuclideanClusterExtraction()
{
}

//----------------------------------------------------------------------------
int vtkPCLEuclideanClusterExtraction::RequestData(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get input and output data objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // perform euclidean cluster extraction
  std::vector<pcl::PointIndices> clusterIndices;
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = vtkPCLConversions::PointCloudFromPolyData(input);
  ApplyEuclideanClusterExtraction(cloud, this->ClusterTolerance, this->MinClusterSize, this->MaxClusterSize, clusterIndices);

  // pass thru input add labels
  vtkSmartPointer<vtkIntArray> labels = vtkPCLConversions::NewLabelsArray(clusterIndices, input->GetNumberOfPoints());
  labels->SetName("cluster_labels");
  output->ShallowCopy(input);
  output->GetPointData()->AddArray(labels);

  return 1;
}

//----------------------------------------------------------------------------
void vtkPCLEuclideanClusterExtraction::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
