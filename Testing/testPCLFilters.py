from vtk import *
from vtkPCLFiltersPython import *


def getSphere(center, radius, resolution=32):
    s = vtkSphereSource()
    s.SetCenter(center)
    s.SetRadius(radius)
    s.SetThetaResolution(resolution)
    s.SetPhiResolution(resolution)
    s.Update()
    return s.GetOutput()


def getPlane(origin, normal, resolution=32):
    p = vtk.vtkPlaneSource()
    p.SetResolution(resolution, resolution)
    p.SetNormal(normal)
    p.SetCenter(origin)
    p.Update()
    return p.GetOutput()


def appendPolyData(*polyData):

    append = vtkAppendPolyData()
    for p in polyData:
        append.AddInput(p)
    append.Update()
    return append.GetOutput()


def getArray(polyData, arrayName):
    return polyData.GetPointData().GetArray(arrayName)


# Test PCL Euclidean Cluster Extraction
s1 = getSphere(center=[-1,0,0], radius=0.5)
s2 = getSphere(center=[1,0,0], radius=0.5)
points = appendPolyData(s1, s2)

f = vtkPCLEuclideanClusterExtraction()
f.SetMinClusterSize(1)
f.SetMaxClusterSize(10000)
f.SetClusterTolerance(0.9)
f.SetInput(points)


def getRange():
    f.Update()
    return f.GetOutput().GetPointData().GetArray('cluster_labels').GetRange()


# Given the center and radii of the spheres, there should be two clusters
assert(getRange() == (1.0, 2.0))

# With this cluster tolerance, there should be only one cluster
f.SetClusterTolerance(1.1)
assert(getRange() == (1.0, 1.0))

# With a very large min size, there should not be any clusters
f.SetMaxClusterSize(int(1e6))
f.SetMinClusterSize(int(1e5))
assert(getRange() == (0.0, 0.0))

# With a very small max size, there should not be any clusters
f.SetMaxClusterSize(5)
f.SetMinClusterSize(1)
assert(getRange() == (0.0, 0.0))


# Test SAC Segmentation Plane filter
plane = getPlane(normal=[1,0,0], origin=[10,0,0])

f = vtkPCLSACSegmentationPlane()
f.SetInput(plane)
f.Update()

print 'plane origin', f.GetPlaneOrigin()
print 'plane normal', f.GetPlaneNormal()
