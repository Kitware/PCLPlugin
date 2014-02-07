import numpy
import vtkPCLFiltersPython as pcl

# perform plane segmentation
f = pcl.vtkPCLSACSegmentationPlane()
f.SetInput(inputs[0].VTKObject)
f.SetDistanceThreshold(0.01)
f.Update()
origin = f.GetPlaneOrigin()
normal = f.GetPlaneNormal()

# for each point, compute signed distance to plane
dist = numpy.dot(inputs[0].Points - origin, normal)

# flip the sign if needed (dot normal with Y axis)
if numpy.dot(normal, [0,1,0]) > 0:
    dist *= -1.0

# pass thru the input data and append the dist_to_plane array
output.ShallowCopy(f.GetOutput())
output.PointData.append(dist, 'dist_to_plane')
