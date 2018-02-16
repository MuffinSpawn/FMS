# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 10:46:59 2018

@author: fms-local
"""
import csv
import math
import numpy
import numpy.linalg as linalg
import matplotlib.pyplot as plt

# Measured Spherical LTCS Coordinates
filename = 'C:\\Users\\fms-local\\Desktop\\FMS\\FMS\\reference_network_LTCS.csv'
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
print('Y (Measured 2D Cartesian LTCS Coordinates):\n{}'.format(Y))
H = numpy.hstack([P[0,2], P[2,2], P[3,2]])-P[0,2]
print('H: {}'.format(H))

# Configured Cylindrical DSCS Coordinates
filename = 'C:\\Users\\fms-local\\Desktop\\FMS\\FMS\\reference_network.csv'
with open(filename, 'r') as file:
    string_data = list(csv.reader(file, delimiter=','))

datap = numpy.ndarray((len(string_data), 3))
for index,row in enumerate(string_data):
    point = list(map(lambda x: float(x), row[1:]))
    datap[index,:] = point
print('Configured Cylindrical DSCS Coordinates:\n{}'.format(datap))

# Configured Cartesian DSCS Coordinates
Pp = numpy.ndarray((5, 3))
for index,point in enumerate(datap):
    Pp[index,0] = point[0] * math.cos(point[1])
    Pp[index,1] = point[0] * math.sin(point[1])
    Pp[index,2] = point[2]
X_all = numpy.vstack([[1, 1, 1, 1, 1], Pp.transpose()])
print('Configured Cartesian DSCS Coordinates:\n{}'.format(X_all))

'''
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

T_inv = linalg.inv(T)
S = numpy.array([1, 0, 0])
Sp_2D = numpy.dot(T_inv, S)
Sp = numpy.array([Sp_2D[1], -dH, Sp_2D[2]])
print('Tracker Cartesian DSCS Coordinates: {}'.format(Sp))

Y = numpy.vstack([[1, 1, 1, 1], numpy.vstack([P[0], P[2], P[4], [0, 0, 0]]).transpose()])
X = numpy.vstack([[1, 1, 1, 1], numpy.vstack([Pp[0], Pp[2], Pp[4], Sp]).transpose()])
T = numpy.dot(Y, linalg.pinv(X))
Y_all = numpy.dot(T, X_all)
print('All Configured LTCS Coordinates (3D LLS):\n{}'.format(Y_all))
'''

A = P[0]
O = numpy.array([0, 0, 0])
AB = P[2]-P[0]
AC = P[4]-P[0]
AO = -P[0]
z = numpy.cross(AB, AC)
z_hat = z / linalg.norm(z)
a = numpy.dot(z_hat, AO)
b = numpy.dot(AB, AO)
c = numpy.dot(AC, AO)

Ap = Pp[0]
ABp = Pp[2]-Pp[0]
ACp = Pp[4]-Pp[0]
zp = numpy.cross(ABp, ACp)
zp_hat = zp / linalg.norm(zp)
M = numpy.vstack([zp_hat, ABp, ACp])
M_inv = linalg.inv(M)
AOp = numpy.dot(M_inv, numpy.array([a, b, c]))
Op = AOp + Ap

Y = numpy.vstack([[1, 1, 1, 1], numpy.vstack([P[0], P[2], P[4], O]).transpose()])
X = numpy.vstack([[1, 1, 1, 1], numpy.vstack([Pp[0], Pp[2], Pp[4], Op]).transpose()])
T = numpy.dot(Y, linalg.pinv(X))
Y_all = numpy.dot(T, X_all)

print('z_hat: {}'.format(z_hat))
print('zp_hat: {}'.format(zp_hat))
print('AB: {}'.format(AB))
print("AB': {}".format(ABp))
print('AC: {}'.format(AC))
print("AC': {}".format(ACp))
print("O': {}".format(Op))
print('All Measured Cartesian LTCS Coordinates:\n{}'.format(P.transpose()))
print('All Configured LTCS Coordinates (3D LLS):\n{}'.format(Y_all))


'''
plt.plot(P[:,2]-P[0,2])
plt.plot(Pp[:,1]-Pp[0,1])
plt.plot(Pp[:,1]-Pp[0,1]+dH)
'''
'''
heights = P[:,2]-numpy.mean(P[:,2])
heightsp = Pp[:,1]-numpy.mean(Pp[:,1])
height_residuals = heights-heightsp
distances = data[:,2]
plt.scatter(distances, height_residuals)
'''
