/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPCLNormalEstimation.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPCLNormalEstimation -
// .SECTION Description
//

#ifndef __vtkPCLNormalEstimation_h
#define __vtkPCLNormalEstimation_h

#include <vtkPolyDataAlgorithm.h>


class vtkPCLNormalEstimation : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkPCLNormalEstimation, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkPCLNormalEstimation *New();

  vtkSetMacro(SearchRadius, double);
  vtkGetMacro(SearchRadius, double);


protected:

  double SearchRadius;

  virtual int RequestData(vtkInformation *request,
                          vtkInformationVector **inputVector,
                          vtkInformationVector *outputVector);

  vtkPCLNormalEstimation();
  virtual ~vtkPCLNormalEstimation();


  virtual int FillInputPortInformation(int port, vtkInformation* info);

private:
  vtkPCLNormalEstimation(const vtkPCLNormalEstimation&);  // Not implemented.
  void operator=(const vtkPCLNormalEstimation&);  // Not implemented.
};

#endif
