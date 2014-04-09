/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPCLNDTRegistration.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPCLNDTRegistration -
// .SECTION Description
//

#ifndef __vtkPCLNDTRegistration_h
#define __vtkPCLNDTRegistration_h

#include <vtkPolyDataAlgorithm.h>
#include <vtkPCLFiltersModule.h>


class VTKPCLFILTERS_EXPORT vtkPCLNDTRegistration : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkPCLNDTRegistration, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkPCLNDTRegistration *New();

  vtkSetMacro(StepSize, double);
  vtkGetMacro(StepSize, double);

  vtkSetMacro(Resolution, double);
  vtkGetMacro(Resolution, double);

  vtkSetMacro(MaxIteration, int);
  vtkGetMacro(MaxIteration, int);

  vtkSetVector3Macro(InitTranslation, double);
  vtkGetVector3Macro(InitTranslation, double);

protected:

  double StepSize;
  double Resolution;
  int MaxIteration;
  double InitTranslation[3];

  virtual int RequestData(vtkInformation *request,
                          vtkInformationVector **inputVector,
                          vtkInformationVector *outputVector);

  vtkPCLNDTRegistration();
  virtual ~vtkPCLNDTRegistration();

private:
  vtkPCLNDTRegistration(const vtkPCLNDTRegistration&);  // Not implemented.
  void operator=(const vtkPCLNDTRegistration&);  // Not implemented.
};

#endif


