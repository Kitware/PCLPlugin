/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPCDReader.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPCDReader.h"

#include "vtkPolyData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkSmartPointer.h"
#include "vtkNew.h"

#include "vtkPCLConversions.h"

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPCDReader);

//----------------------------------------------------------------------------
vtkPCDReader::vtkPCDReader()
{
  this->FileName = 0;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkPCDReader::~vtkPCDReader()
{
  this->SetFileName(NULL);
}

//----------------------------------------------------------------------------
int vtkPCDReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkDataSet *output = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (!this->GetFileName())
    {
    vtkErrorMacro("Filename is not set");
    return 0;
    }

  vtkSmartPointer<vtkPolyData> polyData = vtkPCLConversions::PolyDataFromPCDFile(this->GetFileName());

  if (!polyData)
    {
    vtkErrorMacro("Failed to read pcd file: " << this->GetFileName());
    }

  output->ShallowCopy(polyData);
  return 1;
}

//----------------------------------------------------------------------------
void vtkPCDReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
