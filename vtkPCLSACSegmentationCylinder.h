/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPCLSACSegmentationCylinder.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPCLSACSegmentationCylinder -
// .SECTION Description
//

#ifndef __vtkPCLSACSegmentationCylinder_h
#define __vtkPCLSACSegmentationCylinder_h

#include <vtkPolyDataAlgorithm.h>


class vtkPCLSACSegmentationCylinder : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkPCLSACSegmentationCylinder, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkPCLSACSegmentationCylinder *New();

  vtkSetMacro(DistanceThreshold, double);
  vtkGetMacro(DistanceThreshold, double);

  vtkSetMacro(MaxIterations, int);
  vtkGetMacro(MaxIterations, int);

  vtkSetMacro(RadiusLimit, double);
  vtkGetMacro(RadiusLimit, double);

  vtkSetMacro(SearchRadius, double);
  vtkGetMacro(SearchRadius, double);

  vtkSetMacro(NormalDistanceWeight, double);
  vtkGetMacro(NormalDistanceWeight, double);

protected:

  double NormalDistanceWeight;
  double DistanceThreshold;
  double RadiusLimit;
  double SearchRadius;
  
  int MaxIterations;

  virtual int RequestData(vtkInformation *request,
                          vtkInformationVector **inputVector,
                          vtkInformationVector *outputVector);
  vtkPCLSACSegmentationCylinder();
  virtual ~vtkPCLSACSegmentationCylinder();

private:
  vtkPCLSACSegmentationCylinder(const vtkPCLSACSegmentationCylinder&);  // Not implemented.
  void operator=(const vtkPCLSACSegmentationCylinder&);  // Not implemented.
};

#endif


