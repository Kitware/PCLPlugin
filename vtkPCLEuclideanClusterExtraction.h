/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPCLEuclideanClusterExtraction.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPCLEuclideanClusterExtraction -
// .SECTION Description
//

#ifndef __vtkPCLEuclideanClusterExtraction_h
#define __vtkPCLEuclideanClusterExtraction_h

#include <vtkPolyDataAlgorithm.h>


class vtkPCLEuclideanClusterExtraction : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkPCLEuclideanClusterExtraction, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkPCLEuclideanClusterExtraction *New();

  vtkSetMacro(ClusterTolerance, double);
  vtkGetMacro(ClusterTolerance, double);

  vtkSetMacro(MinClusterSize, int);
  vtkGetMacro(MinClusterSize, int);

  vtkSetMacro(MaxClusterSize, int);
  vtkGetMacro(MaxClusterSize, int);

protected:

  double ClusterTolerance;
  int MinClusterSize;
  int MaxClusterSize;

  virtual int RequestData(vtkInformation *request,
                          vtkInformationVector **inputVector,
                          vtkInformationVector *outputVector);

  vtkPCLEuclideanClusterExtraction();
  virtual ~vtkPCLEuclideanClusterExtraction();

private:
  vtkPCLEuclideanClusterExtraction(const vtkPCLEuclideanClusterExtraction&);  // Not implemented.
  void operator=(const vtkPCLEuclideanClusterExtraction&);  // Not implemented.
};

#endif


