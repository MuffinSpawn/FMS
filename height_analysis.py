# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 10:46:59 2018

@author: Peter G. Lane (petergwinlane@gmail.com)
"""
import csv
import math
import numpy
import numpy.linalg as linalg
# import matplotlib.pyplot as plt

def spherical_to_cartesian(spherical_coordinates):
    points = spherical_coordinates.transpose()
    cartesian_coordinates = numpy.ndarray(numpy.shape(points))
    cartesian_coordinates[0,:] = points[0] * numpy.cos(points[1]) * numpy.sin(points[2])
    cartesian_coordinates[1,:] = points[0] * numpy.sin(points[1]) * numpy.sin(points[2])
    cartesian_coordinates[2,:] = points[0] * numpy.cos(points[2])
    return cartesian_coordinates.transpose()

def cartesian_to_spherical(cartesian_coordinates):
    points = cartesian_coordinates.transpose()
    cylindrical_coordinates = numpy.ndarray(numpy.shape(points))
    cylindrical_coordinates[0] = numpy.sqrt(numpy.sum(points**2), axis=1)
    cylindrical_coordinates[1] = numpy.arctan(points[1], points[0])
    cylindrical_coordinates[2] = numpy.arccos(points[2]/cylindrical_coordinates[0])
    return cylindrical_coordinates.transpose()

def cylindrical_to_cartesian(cylindrical_coordinates):
    points = cylindrical_coordinates.transpose()
    cartesian_coordinates = numpy.ndarray(numpy.shape(points))
    cartesian_coordinates[0] = points[0] * numpy.cos(points[1])
    cartesian_coordinates[1] = points[0] * numpy.sin(points[1])
    cartesian_coordinates[2] = points[2]
    return cartesian_coordinates.transpose()

def cartesian_to_cylindrical(cartesian_coordinates):
    points = cartesian_coordinates.transpose()
    cylindrical_coordinates = numpy.ndarray(numpy.shape(points))
    cylindrical_coordinates[0] = numpy.sqrt(points[0]**2 + points[1]**2)
    cylindrical_coordinates[1] = numpy.arctan2(points[1], points[0])
    cylindrical_coordinates[2] = points[2]
    return cylindrical_coordinates.transpose()

def calculate_transform(spherical_LTCS_coordinates, cylindrical_DSCS_coordinates):
    cartesian_LTCS_coordinates = spherical_to_cartesian(spherical_LTCS_coordinates)
    cartesian_DSCS_coordinates = cylindrical_to_cartesian(cylindrical_DSCS_coordinates)

    # calculate the LTCS-to-DSCS transform matrix
    constant_vector = numpy.ones(numpy.shape(cartesian_LTCS_coordinates)[0])
    X = numpy.vstack([constant_vector, cartesian_LTCS_coordinates.transpose()])
    Y = numpy.vstack([constant_vector, cartesian_DSCS_coordinates.transpose()])
    print('X:\n{}'.format(X))
    print('Y:\n{}'.format(Y))
    print('pinv(x):{}'.format(linalg.pinv(X)))
    transform = numpy.dot(Y, linalg.pinv(X))
    print('Transform: {}'.format(transform))
    
    projection_matrix = numpy.dot(linalg.pinv(X), X)
    print('Projection ("Hat") Matrix:\n{}'.format(projection_matrix))
    

    Yp = numpy.dot(transform, X)
    print('Yp: {}'.format(Yp))
    
    means = numpy.mean(Y, axis=1)
    print('Mean: {}'.format(means))
    total_sum_of_squares = numpy.sum((Y.transpose() - means)**2, axis=0)
    total_sum_of_squares[0] = 1.0e12
    print('TSS: {}'.format(total_sum_of_squares))
    
    residuals = Y - Yp
    print('Residuals: {}'.format(residuals))
    ''' The normal R^2 is not expressive enough to indicate bad mappings.
    residual_sum_of_squares = numpy.sum(residuals**2, axis=1)
    print('RSS: {}'.format(residual_sum_of_squares))
    # r_squared = 1 - residual_sum_of_squares / total_sum_of_squares
    '''

    print('1-leverage: {}'.format(1-projection_matrix.diagonal()))
    predictive_residuals = residuals/(1-projection_matrix.diagonal())
    print('Pred. Res.: {}'.format(predictive_residuals))
    predictive_residual_sum_of_squares = numpy.sum(predictive_residuals**2, axis=1)
    print('PRESS: {}'.format(predictive_residual_sum_of_squares))
    predictive_r_squared = \
        1 - predictive_residual_sum_of_squares/total_sum_of_squares

    return (transform, predictive_r_squared)

def LTCS_to_DSCS(transform, spherical_LTCS_coordinates):
    cartesian_LTCS_coordinates = spherical_to_cartesian(spherical_LTCS_coordinates)
    constant_vector = numpy.ones(numpy.shape(cartesian_LTCS_coordinates)[0])
    X = numpy.vstack([constant_vector, cartesian_LTCS_coordinates.transpose()])
    Y = numpy.dot(transform, X)
    return(cartesian_to_cylindrical(Y[1:].transpose()))

def DSCS_to_LTCS(transform, cylindrical_LTCS_coordinates):
    cartesian_DSCS_coordinates = cylindrical_to_cartesian(cylindrical_LTCS_coordinates)
    constant_vector = numpy.ones(numpy.shape(cartesian_DSCS_coordinates)[0])
    Y = numpy.vstack([constant_vector, cartesian_DSCS_coordinates.transpose()])
    X = numpy.dot(linalg.inv(transform), Y)
    return(cartesian_to_spherical(X[1:].transpose()))

# Measured Spherical LTCS Coordinates
#filename = 'reference_network_LTCS.csv'
#filename = 'reference_network_LTCS_20180504.csv'
#filename = 'reference_network_LTCS_TEST.csv'
filename = 'reference_network_IB2_LTCS.csv'
with open(filename, 'r') as file:
    string_data = list(csv.reader(file, delimiter=','))

data_LTCS = numpy.ndarray((len(string_data), 3))
for index,row in enumerate(string_data):
    point = list(map(lambda x: float(x), row[1:]))
    data_LTCS[index,:] = point
print('Measured Spherical LTCS Coordinates (azimuth, zenith, distance):\n{}'\
      .format(data_LTCS))

# Measured Cartesian LTCS Coordinates
P = spherical_to_cartesian(data_LTCS)
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
filename = 'reference_network_IB2.csv'
with open(filename, 'r') as file:
    string_data = list(csv.reader(file, delimiter=','))

data_DSCS = numpy.ndarray((len(string_data), 3))
for index,row in enumerate(string_data):
    point = list(map(lambda x: float(x), row[1:]))
    data_DSCS[index,:] = point
print('Configured Cylindrical DSCS Coordinates (rho, theta, z):\n{}'.format(data_DSCS))

# Configured Cartesian DSCS Coordinates
Pp = cylindrical_to_cartesian(data_DSCS)
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

transform,r_squared = calculate_transform(data_LTCS, data_DSCS)
print('R^2 of LTCS->DSCS Transform: {}'.format(r_squared))

cartesian_data_DSCS = cylindrical_to_cartesian(data_DSCS)
transformed_DSCS = LTCS_to_DSCS(transform, data_LTCS)
cartesian_transformed_DSCS = cylindrical_to_cartesian(transformed_DSCS)
print(cartesian_data_DSCS)
print(cartesian_transformed_DSCS)
cartesian_residuals_DSCS = cartesian_data_DSCS - cartesian_transformed_DSCS
print('Cartesian DSCS Residuals: {}'.format(cartesian_residuals_DSCS))
