# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 10:46:59 2018

@author: fms-local
"""
import csv
import math
import numpy
import numpy.linalg as linalg

filename = 'C:\\Users\\fms-local\\Desktop\\FMS\\reference_network.csv'
with open(filename, 'r') as file:
    string_data = list(csv.reader(file, delimiter=','))

data = numpy.ndarray((len(string_data), 3))
for index,row in enumerate(string_data):
    point = list(map(lambda x: float(x), row[1:]))
    data[index,:] = point
print('Config Coordinates: {}\n'.format(data))


P1_sph = numpy.array([2481.968, 0.050865, 1.51838])
P2_sph = numpy.array([4162.099, 2.351501, 1.56747])
P3_sph = numpy.array([1.011E+4, 2.532285, 1.57984])

P1 = numpy.array([P1_sph[0] * math.cos(P1_sph[1]) * math.sin(P1_sph[2]),\
                 P1_sph[0] * math.sin(P1_sph[1]) * math.sin(P1_sph[2]),\
                 P1_sph[0] * math.cos(P1_sph[2])])
P2 = numpy.array([P2_sph[0] * math.cos(P2_sph[1]) * math.sin(P2_sph[2]),\
                 P2_sph[0] * math.sin(P2_sph[1]) * math.sin(P2_sph[2]),\
                 P2_sph[0] * math.cos(P2_sph[2])])
P3 = numpy.array([P3_sph[0] * math.cos(P3_sph[1]) * math.sin(P3_sph[2]),\
                 P3_sph[0] * math.sin(P3_sph[1]) * math.sin(P3_sph[2]),\
                 P3_sph[0] * math.cos(P3_sph[2])])
P0 = numpy.array([0., 0., 0.])
X = numpy.vstack([P1, P2, P3, P0]).transpose()
print('X (Cartesian LTCS):\n{}'.format(X))

# (2)
P12 = P2 - P1
P13 = P3 - P1
P10 = P0 - P1

# (3)
P10_norm = linalg.norm(P1)
print("|P10|: {}".format(P10_norm))

# (4)
z = numpy.cross(P12, P13)
z_hat = z / linalg.norm(z)
print('z_hat: {}'.format(z_hat))
print(numpy.dot(P12, z_hat))
print(numpy.dot(P13, z_hat))

# (5)
a = numpy.dot(z_hat, P10)
print('a: {}'.format(a))
b = numpy.dot(P12, P10)
c = numpy.dot(P13, P10)

P1p = numpy.array([data[0,0] * math.cos(data[0,1]),\
                  data[0,0] * math.sin(data[0,1]),\
                  data[0,2]])
P2p = numpy.array([data[1,0] * math.cos(data[1,1]),\
                  data[1,0] * math.sin(data[1,1]),\
                  data[1,2]])
P3p = numpy.array([data[2,0] * math.cos(data[2,1]),\
                  data[2,0] * math.sin(data[2,1]),\
                  data[2,2]])

P12p = P2p - P1p
P13p = P3p - P1p
zp = numpy.cross(P12p, P13p)
z_hatp = zp / linalg.norm(zp)
print("z_hat': {}".format(z_hatp))

M = numpy.vstack([z_hatp, P12p, P13p])
print('M:\n{}'.format(M))
M_inv = linalg.inv(M)
print('M_inv:\n{}'.format(M_inv))
P12c = numpy.array([a, b, c])
print('P12c: {}'.format(P12c))
P1P0p = numpy.dot(M_inv, P12c)
print("P10': {}".format(P1P0p))
print("|P10'|: {}".format(linalg.norm(P1P0p)))
P0p = P1P0p + P1p
print("O': {}".format(P0p))

Y = numpy.vstack([P1p, P2p, P3p, P0p]).transpose()
print('Y (Cartesian DSCS):\n{}'.format(Y))

# LTCS (Spherical) Coordinates
#
# Name, theta, phi, r
# B2, 0.050865, 1.51838, 2481.968
# B3, 2.351501, 1.56747, 4162.099
# B4, 2.532285, 1.57984, 1.011E+4
# B6, -2.617606, 1.55448, 1.868E+4
# B8, -1.918747, 1.55075, 1.593E+4
# 1, -2.136088, 1.79359, 3463.905