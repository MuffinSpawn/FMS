# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 10:46:59 2018

@author: fms-local
"""
import csv
import math
import numpy
import numpy.linalg as linalg

# Measured Spherical LTCS Coordinates
filename = 'C:\\Users\\fms-local\\Desktop\\FMS\\reference_network_LTCS.csv'
with open(filename, 'r') as file:
    string_data = list(csv.reader(file, delimiter=','))

data = numpy.ndarray((len(string_data), 3))
for index,row in enumerate(string_data):
    point = list(map(lambda x: float(x), row[1:]))
    data[index,:] = point
print('Measured Spherical LTCS Coordinates:\n{}'.format(data))

# Measured Cartesian LTCS Coordinates
P = numpy.ndarray((5, 3))
for index,point in enumerate(data):
    P[index,0] = point[2] * math.cos(point[0]) * math.sin(point[1])
    P[index,1] = point[2] * math.sin(point[0]) * math.sin(point[1])
    P[index,2] = point[2] * math.cos(point[1])
print('Measured Cartesian LTCS Coordinates:\n{}'.format(P.transpose()))

Y = numpy.vstack([[1, 1, 1], numpy.vstack([P[0,:2], P[2,:2], P[4,:2]]).transpose()])
print('Y (Measured 2D Cartesian LTCS Coordinates):\n{}'.format(X))
H = numpy.hstack([P[0,2], P[2,2], P[3,2]])-P[0,2]
print('H: {}'.format(H))

# Configured Cylindrical DSCS Coordinates
filename = 'C:\\Users\\fms-local\\Desktop\\FMS\\reference_network.csv'
with open(filename, 'r') as file:
    string_data = list(csv.reader(file, delimiter=','))

data = numpy.ndarray((len(string_data), 3))
for index,row in enumerate(string_data):
    point = list(map(lambda x: float(x), row[1:]))
    data[index,:] = point
print('Configured Cylindrical DSCS Coordinates:\n{}'.format(data))

# Configured Cartesian DSCS Coordinates
Pp = numpy.ndarray((5, 3))
for index,point in enumerate(data):
    Pp[index,0] = point[0] * math.cos(point[1])
    Pp[index,1] = point[0] * math.sin(point[1])
    Pp[index,2] = point[2]
X_all = numpy.vstack([[1, 1, 1, 1, 1], Pp.transpose()])
print('Configured Cartesian DSCS Coordinates:\n{}'.format(X_all))

mask = [True, False, True]
X = numpy.vstack([[1, 1, 1], numpy.vstack([Pp[0,mask], Pp[2,mask], Pp[4,mask]]).transpose()])
print('X (Configured 2D Cartesian LTCS Coordinates):\n{}'.format(Y))
Hp = numpy.hstack([Pp[0,1], Pp[2,1], Pp[3,1]])-Pp[0,1]
print('Hp: {}'.format(Hp))
Hp_all = Pp[:,1]

T = numpy.dot(Y, linalg.pinv(X))
print('DSCS-LTCSS Transform Matrix:\n{}'.format(T))
dH = numpy.mean(H - Hp)
print('Height Residuals: {}'.format(H-Hp))
print('dH: {}'.format(dH))
print(Hp_all+dH)

print('All Measured Cartesian LTCS Coordinates:\n{}'.format(P.transpose()))

Y_all_2D = numpy.dot(T, X_all[(True, True, False, True),:])
Y_all = numpy.vstack([Y_all_2D[1], Y_all_2D[2], Hp_all+dH])
print('All Configured LTCS Coordinates:\n{}'.format(Y_all))

''' 3D LLS doesn't work with only 3 points
Y = numpy.vstack([[1, 1, 1], numpy.vstack([P[0], P[2], P[4]]).transpose()])
X = numpy.vstack([[1, 1, 1], numpy.vstack([Pp[0], Pp[2], Pp[4]]).transpose()])
T = numpy.dot(Y, linalg.pinv(X))
Y_all = numpy.dot(T, X_all)
print('All Configured LTCS Coordinates (3D LLS):\n{}'.format(Y_all))
'''
