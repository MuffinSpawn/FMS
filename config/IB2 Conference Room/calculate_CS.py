# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 12:17:43 2018

@author: plane
"""

import math
import pandas
import numpy
import numpy.linalg as linalg

def spherical_to_cartesian(coordinates_SCC):
    coordinates_RHR = numpy.ndarray(numpy.shape(coordinates_SCC), dtype=numpy.float64)
    for index,point in enumerate(coordinates_SCC):
        coordinates_RHR[index,0] = point[0]*math.cos(point[1])*math.sin(point[2])
        coordinates_RHR[index,1] = point[0]*math.sin(point[1])*math.sin(point[2])
        coordinates_RHR[index,2] = point[0]*math.cos(point[2])
    return coordinates_RHR

def cartesian_to_cylindrical(cartesian_DSCS):
    cylindrical_DSCS = numpy.ndarray(numpy.shape(cartesian_DSCS))
    for index,cartesian_point in enumerate(cartesian_DSCS):
        cylindrical_DSCS[index,0] = math.sqrt(numpy.sum(cartesian_point[:2]**2))
        cylindrical_DSCS[index,1] = math.atan2(cartesian_point[1], cartesian_point[0])
        cylindrical_DSCS[index,2] = cartesian_point[2]
    return cylindrical_DSCS


data_spherical = pandas.read_csv('mapping_LTCS.csv', dtype=numpy.float64)
data_cartesian = spherical_to_cartesian(data_spherical.as_matrix())

# data['Theta [rad]'] = list(map(lambda x: x/1000-2*math.pi, data['Theta [mrad]']))
# data.to_csv('argonne_ref_network.csv', index=False, float_format='%.6g', header=False, columns=('Point #', 'R [mm]', 'Theta [rad]', 'Z[mm]'))

centroid = numpy.mean(data_cartesian, axis=0)
data_translated = data_cartesian - centroid
U,s,v = linalg.svd(data_translated)

z = v[0]
z[2] = 0.0
z_hat_LTCS = z / linalg.norm(z)
origin_LTCS = centroid - 1500*z_hat_LTCS
y_hat_LTCS = numpy.array([0.0, 0.0, 1.0])
x_hat_LTCS = numpy.cross(y_hat_LTCS, z_hat_LTCS)

points_LTCS = numpy.vstack((origin_LTCS, origin_LTCS+1000*x_hat_LTCS, origin_LTCS+1000*y_hat_LTCS, origin_LTCS+1000*z_hat_LTCS))
points_DSCS = numpy.array([[0.0,0.0,0.0], [1000.0,0.0,0.0], [0.0,1000.0,0.0],[0.0,0.0,1000.0]])

'''
    # calculate the LTCS-to-DSCS transform matrix
    X = numpy.vstack([[1, 1, 1, 1], numpy.vstack([A, B, C, S]).transpose()])
    Y = numpy.vstack([[1, 1, 1, 1], numpy.vstack([Ap, Bp, Cp, Sp]).transpose()])
    return numpy.dot(Y, linalg.pinv(X))
'''
X = numpy.vstack(([1,1,1,1], points_LTCS.transpose()))
Y = numpy.vstack(([1,1,1,1], points_DSCS.transpose()))
transform_LTCS_DSCS = numpy.dot(Y, linalg.pinv(X))

ref_spherical = pandas.read_csv('reference_network_LTCS.csv', dtype=numpy.float64)
ref_cartesian_LTCS = spherical_to_cartesian(ref_spherical.as_matrix())
X_ref = numpy.vstack(([1,1,1], ref_cartesian_LTCS.transpose()))
ref_cartesian_DSCS = numpy.dot(transform_LTCS_DSCS, X_ref)
ref_cylindrical_DSCS = cartesian_to_cylindrical(ref_cartesian_DSCS[1:].transpose())
print(ref_cylindrical_DSCS)

origin_DSCS = numpy.array([[1.0, 0.0, 0.0, 0.0]]).transpose()
print(numpy.dot(transform_LTCS_DSCS, origin_DSCS))