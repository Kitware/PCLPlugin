/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkBoundsExtentTranslator.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkBoundsExtentTranslator - Extent translation through lookup table.
// .SECTION Description
// vtkBoundsExtentTranslator provides a vtkExtentTranslator that is
// programmed with a specific extent corresponding to each piece
// number.  Readers can provide this to an application to allow the
// pipeline to execute using the same piece breakdown that is provided
// in the input file.

#ifndef __vtkBoundsExtentTranslator_h
#define __vtkBoundsExtentTranslator_h

#include "vtkExtentTranslator.h"
#include <vector>

class VTK_EXPORT vtkBoundsExtentTranslator : public vtkExtentTranslator
{
public:
  vtkTypeMacro(vtkBoundsExtentTranslator,vtkExtentTranslator);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  static vtkBoundsExtentTranslator* New();

  // Description:
  // Set the number of pieces that will be stored.
  // Executives set this and the value must match that set by the
  // PartitionFilter that generated this data, otherwise failure
  // is probable
  virtual void SetNumberOfPieces(int pieces);

  // Description:
  // Used by executives to compute a structured extent.
  // Only valid when a default Spacing has been set, otherwise
  // all pieces return the WholeExtent
  virtual int PieceToExtent();  
  virtual int PieceToExtentByPoints();
  virtual int BoundsToExtentThreadSafe(
      double *bounds, int *wholeExtent, int *resultExtent);
  virtual int PieceToExtentThreadSafe(int piece, int numPieces, 
                              int ghostLevel, int *wholeExtent, 
                              int *resultExtent, int splitMode, 
                              int byPoints);
  
  // Description:
  // Set the extent to be used for a piece.  This sets the extent table
  // entry for the piece.
  virtual void SetBoundsForPiece(int piece, double* bounds);
  
  // Description:  
  // Get the bounds table entry for the given piece.  
  // Structured Extents should be calculated using PieceToExtent
  // (after a valid spacing has been set).
  virtual void    GetBoundsForPiece(int piece, double *bounds);
  virtual double *GetBoundsForPiece(int piece);
  
  // Description:
  // Set the maximum ghost overlap region that is required 
  vtkSetMacro(MaximumGhostDistance, double);
  vtkGetMacro(MaximumGhostDistance, double);
  
  // Description:
  // To convert a bounding box to a StructuredExtent, a spacing
  // (between samples) is necessary - and the global bounds.
  vtkSetVector3Macro(Spacing, double);
  vtkGetVector3Macro(Spacing, double);

  // Description:
  // After setting all bounds, call this to compute the global/whole bounds
  virtual void InitWholeBounds();

  // Description:
  // Returns the union of all bounds
  virtual double *GetWholeBounds();

protected:
   vtkBoundsExtentTranslator();
  ~vtkBoundsExtentTranslator();
  
  // Store the extent table in a single array.  Every 6 values form an extent.
  std::vector<double> BoundsTable;
  double WholeBounds[6];
  double Spacing[3];
  double MaximumGhostDistance;
   
private:
  vtkBoundsExtentTranslator(const vtkBoundsExtentTranslator&);  // Not implemented.
  void operator=(const vtkBoundsExtentTranslator&);  // Not implemented.
};

#endif
