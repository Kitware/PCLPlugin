/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPCDReader.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPCDReader -
// .SECTION Description
//

#ifndef __vtkPCDReader_h
#define __vtkPCDReader_h

#include <vtkPolyDataAlgorithm.h>
#include <vtkPCLFiltersModule.h>

class VTKPCLFILTERS_EXPORT vtkPCDReader : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkPCDReader, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkPCDReader *New();

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

protected:

  char* FileName;

  virtual int RequestData(vtkInformation *request,
                          vtkInformationVector **inputVector,
                          vtkInformationVector *outputVector);


  vtkPCDReader();
  virtual ~vtkPCDReader();

private:
  vtkPCDReader(const vtkPCDReader&);  // Not implemented.
  void operator=(const vtkPCDReader&);  // Not implemented.
};

#endif


