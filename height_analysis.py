# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 10:46:59 2018

@author: fms-local
"""
import csv
import math
import numpy
import numpy.linalg as linalg
# import matplotlib.pyplot as plt

# Measured Spherical LTCS Coordinates
filename = 'C:\\Users\\fms-local\\Desktop\\FMS\\reference_network_LTCS.csv'
with open(filename, 'r') as file:
    string_data = list(csv.reader(file, delimiter=','))

data = numpy.ndarray((len(string_data), 3))
for index,row in enumerate(string_data):
    point = list(map(lambda x: float(x), row[1:]))
    data[index,:] = point
print('Measured Spherical LTCS Coordinates (azimuth, zenith, distance):\n{}'\
      .format(data.transpose()))

# Measured Cartesian LTCS Coordinates
P = numpy.ndarray((5, 3))
for index,point in enumerate(data):
    P[index,0] = point[2] * math.cos(point[0]) * math.sin(point[1])
    P[index,1] = point[2] * math.sin(point[0]) * math.sin(point[1])
    P[index,2] = point[2] * math.cos(point[1])
print('Measured Cartesian LTCS Coordinates (x, y, z (up)):\n{}'.format(P.transpose()))
mean_height = numpy.mean(P[:,2])
# relative_heights = mean_height-P[:,2]
relative_heights = P[:,2]-P[0,2]
print('Measured Relative Heights: {}'.format(relative_heights))
print('Measured Height Mean: {}'.format(numpy.mean(relative_heights)))
print('Measured Height Std. Dev.: {}'.format(numpy.std(relative_heights)))
centroid = numpy.array([numpy.mean(P[:,0]), numpy.mean(P[:,1])])
print('Measured Horizontal Centroid: {}'.format(centroid))
radii = numpy.sqrt((centroid[0]-P[:,0])**2 + (centroid[1]-P[:,1])**2)
print('Measured Horizontal Radii: {}'.format(radii))
print('Measured Horizontal Radius Mean: {}'.format(numpy.mean(radii)))
print('Measured Horizontal Radius Std. Dev.: {}'.format(numpy.std(radii)))
relative_distances = numpy.sqrt((P[:,0]-P[0,0])**2 + (P[:,1]-P[0,1])**2 + (P[:,2]-P[0,2])**2)
print('Measured Relative (B2) Distances: {}'.format(relative_distances))
print()

# Configured Cylindrical DSCS Coordinates
filename = 'C:\\Users\\fms-local\\Desktop\\FMS\\reference_network.csv'
with open(filename, 'r') as file:
    string_data = list(csv.reader(file, delimiter=','))

data = numpy.ndarray((len(string_data), 3))
for index,row in enumerate(string_data):
    point = list(map(lambda x: float(x), row[1:]))
    data[index,:] = point
print('Configured Cylindrical DSCS Coordinates (rho, theta, z):\n{}'.format(data))

# Configured Cartesian DSCS Coordinates
Pp = numpy.ndarray((5, 3))
for index,point in enumerate(data):
    Pp[index,0] = point[0] * math.cos(point[1])
    Pp[index,1] = point[0] * math.sin(point[1])
    Pp[index,2] = point[2]
X_all = Pp.transpose()
print('Configured Cartesian DSCS Coordinates: (x, y (up), z)\n{}'.format(X_all))
mean_height = numpy.mean(Pp[:,1])
# relative_heights = mean_height-Pp[:,1]
relative_heights = Pp[:,1]-Pp[0,1]
print('Configured Relative Heights: {}'.format(relative_heights))
print('Configured Height Mean: {}'.format(numpy.mean(relative_heights)))
print('Configured Height Std. Dev.: {}'.format(numpy.std(relative_heights)))
centroidp = numpy.array([numpy.mean(Pp[:,0]), numpy.mean(Pp[:,2])])
print('Configured Horizontal Centroid: {}'.format(centroidp))
radiip = numpy.sqrt((centroidp[0]-Pp[:,0])**2 + (centroidp[1]-Pp[:,2])**2)
print('Configured Horizontal Radii: {}'.format(radiip))
print('Configured Horizontal Radius Mean: {}'.format(numpy.mean(radiip)))
print('Configured Horizontal Radius Std. Dev.: {}'.format(numpy.std(radiip)))
relative_distances = numpy.sqrt((Pp[:,0]-Pp[0,0])**2 + (Pp[:,1]-Pp[0,1])**2 + (Pp[:,2]-Pp[0,2])**2)
print('Configured Relative (B2) Distances: {}'.format(relative_distances))
print()

height_residuals = Pp[:,1]-P[:,2]
print('Height Residuals: {}'.format(height_residuals))
print('Height Residuals Mean: {}'.format(numpy.mean(height_residuals)))
print('Height Residuals Std. Dev.: {}'.format(numpy.std(height_residuals)))
radii_residuals = radiip-radii
print('Radii Residuals: {}'.format(radii_residuals))
print('Radii Residuals Mean: {}'.format(numpy.mean(radii_residuals)))
print('Radii Residuals Std. Dev.: {}'.format(numpy.std(radii_residuals)))
