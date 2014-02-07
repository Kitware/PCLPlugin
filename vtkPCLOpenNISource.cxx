/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPCLOpenNISource.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPCLOpenNISource.h"
#include "vtkPCLConversions.h"

#include "vtkPolyData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkSmartPointer.h"

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/openni_grabber.h>
#include <boost/thread/thread.hpp>

typedef pcl::PointCloud<pcl::PointXYZRGBA> Cloud;
typedef Cloud::Ptr CloudPtr;
typedef Cloud::ConstPtr CloudConstPtr;

//----------------------------------------------------------------------------
class vtkPCLOpenNISource::vtkInternal
{
public:

  vtkInternal()
  {
    this->Grabber = 0;
    this->NewData = false;
  }

  ~vtkInternal()
  {
    delete this->Grabber; 
  }

  void HandleIncomingCloud(const CloudConstPtr& newCloud)
  {
    vtkSmartPointer<vtkPolyData> newPolyData = vtkPCLConversions::PolyDataFromPointCloud(newCloud);
    boost::lock_guard<boost::mutex> lock(this->mutex);
    this->PolyData = newPolyData;
    this->NewData = true;
  }

  vtkSmartPointer<vtkPolyData> GetLatestPolyData()
  {
    boost::lock_guard<boost::mutex> lock(this->mutex);
    vtkSmartPointer<vtkPolyData> polyData = this->PolyData;
    this->PolyData = NULL;
    this->NewData = false;
    return polyData;
  }

  bool HasNewData()
  {
    boost::lock_guard<boost::mutex> lock(this->mutex);
    return this->NewData;
  }

  bool NewData;
  pcl::OpenNIGrabber* Grabber;
  boost::mutex mutex;
  vtkSmartPointer<vtkPolyData> PolyData;

  boost::function<void (const CloudConstPtr&)> Callback;

};

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPCLOpenNISource);

//----------------------------------------------------------------------------
vtkPCLOpenNISource::vtkPCLOpenNISource()
{
  this->Internal = new vtkInternal;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkPCLOpenNISource::~vtkPCLOpenNISource()
{
  delete this->Internal;
}

//----------------------------------------------------------------------------
void vtkPCLOpenNISource::StartGrabber()
{
  if (!this->Internal->Grabber)
    {
    this->Internal->Grabber = new pcl::OpenNIGrabber("");
    this->Internal->Callback = boost::bind(&vtkPCLOpenNISource::vtkInternal::HandleIncomingCloud, this->Internal, _1);
    this->Internal->Grabber->registerCallback(this->Internal->Callback);
    }
  this->Internal->Grabber->start();
}

//----------------------------------------------------------------------------
void vtkPCLOpenNISource::StopGrabber()
{
  if (this->Internal->Grabber)
    {
    this->Internal->Grabber->stop();
    }
}

//----------------------------------------------------------------------------
bool vtkPCLOpenNISource::HasNewData()
{
  return this->Internal->HasNewData();
}

//----------------------------------------------------------------------------
void vtkPCLOpenNISource::Poll()
{
  if (this->HasNewData())
    {
    this->Modified();
    }
}

//----------------------------------------------------------------------------
int vtkPCLOpenNISource::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkDataSet *output = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (!this->HasNewData())
    {
    return 1;
    }

  output->ShallowCopy(this->Internal->GetLatestPolyData());
  return 1;
}

//----------------------------------------------------------------------------
void vtkPCLOpenNISource::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
