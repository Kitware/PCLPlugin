/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPCLNDTRegistration.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPCLNDTRegistration.h"
#include "vtkPCLConversions.h"

#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkSmartPointer.h"
#include "vtkNew.h"

#include <pcl/registration/ndt.h>
#include <pcl/filters/approximate_voxel_grid.h>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPCLNDTRegistration);

//----------------------------------------------------------------------------
vtkPCLNDTRegistration::vtkPCLNDTRegistration()
{
  this->StepSize = .1;
  this->Resolution = 1.0;
  this->MaxIteration = 35;

  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkPCLNDTRegistration::~vtkPCLNDTRegistration()
{
}

//----------------------------------------------------------------------------
int vtkPCLNDTRegistration::RequestData(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get input and output data objects
  vtkInformation *inInfo1 = inputVector[0]->GetInformationObject(0);
  vtkPolyData *input1 = vtkPolyData::SafeDownCast(inInfo1->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *inInfo2 = inputVector[1]->GetInformationObject(0);
  vtkPolyData *input2 = vtkPolyData::SafeDownCast(inInfo2->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // perform outlier removal
  pcl::PointCloud<pcl::PointXYZ>::Ptr source_cloud = vtkPCLConversions::PointCloudFromPolyData(input1);

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2 = vtkPCLConversions::PointCloudFromPolyData(input2);

  // Filtering input scan to roughly 10% of original size to increase speed of registration.
  pcl::PointCloud<pcl::PointXYZ>::Ptr filtered_cloud (new pcl::PointCloud<pcl::PointXYZ>);
  pcl::ApproximateVoxelGrid<pcl::PointXYZ> approximate_voxel_filter;
  approximate_voxel_filter.setLeafSize (this->Resolution * .2,
                                        this->Resolution * .2,
                                        this->Resolution * .2);
  approximate_voxel_filter.setInputCloud (source_cloud);
  approximate_voxel_filter.filter (*filtered_cloud);
  std::cout << "Filtered cloud contains " << filtered_cloud->size ()
            << " data points from room_scan2.pcd" << std::endl;


  pcl::NormalDistributionsTransform<pcl::PointXYZ, pcl::PointXYZ> ndt;
  ndt.setStepSize(this->StepSize);
  ndt.setResolution(this->Resolution);
  ndt.setMaximumIterations(this->MaxIteration);
  ndt.setTransformationEpsilon (0.01);

  ndt.setInputSource(filtered_cloud);
  ndt.setInputTarget(cloud2);

  Eigen::AngleAxisf init_rotation (0.0, Eigen::Vector3f::UnitZ ());
  Eigen::Translation3f inittranslate(this->InitTranslation[0],
                                        this->InitTranslation[1],
                                        this->InitTranslation[2]);
  Eigen::Matrix4f init = (inittranslate * init_rotation).matrix();

  pcl::PointCloud<pcl::PointXYZ>::Ptr output_cloud(new pcl::PointCloud<pcl::PointXYZ>);
  ndt.align(*output_cloud, init);

  pcl::transformPointCloud(*source_cloud, *output_cloud, ndt.getFinalTransformation());

  // pass thru input add labels
  //vtkSmartPointer
  output->DeepCopy(input1);
  vtkSmartPointer<vtkPolyData> outputPoly = vtkPCLConversions::PolyDataFromPointCloud(output_cloud);
  output->GetPoints()->ShallowCopy(outputPoly->GetPoints());

  return 1;
}

//----------------------------------------------------------------------------
void vtkPCLNDTRegistration::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
