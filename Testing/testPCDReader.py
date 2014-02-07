# This script reads a PCD files and runs the
# PCL to vtkPolyData conversion benchmark

import sys
if len(sys.argv) < 2:
    print 'Usage: %s <filename.pcd>' % sys.argv[0]
    sys.exit(1)

from vtkPCLFiltersPython import *

pcdFile = sys.argv[1]

r = vtkPCDReader()
r.SetFileName(pcdFile)
r.Update()

p = r.GetOutput()
vtkPCLConversions.PerformPointCloudConversionBenchmark(p)
