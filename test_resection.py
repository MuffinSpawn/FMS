# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 10:46:59 2018

@author: fms-local
"""
import csv
import math
import numpy
import numpy.linalg as linalg

# Measured Cylindrical DSCS Coordinates
filename = 'C:\\Users\\fms-local\\Desktop\\FMS\\reference_network.csv'
with open(filename, 'r') as file:
    string_data = list(csv.reader(file, delimiter=','))

data = numpy.ndarray((len(string_data), 3))
for index,row in enumerate(string_data):
    point = list(map(lambda x: float(x), row[1:]))
    data[index,:] = point
print('Config Coordinates: {}\n'.format(data))

# Measured Spherical LTCS Coordinates
P1 = numpy.array([2481.968, 0.050865, 1.51838])
P2 = numpy.array([4162.099, 2.351501, 1.56747])
P3 = numpy.array([1.011E+4, 2.532285, 1.57984])
P0 = numpy.array([0., 0., 0.])
X = numpy.vstack([P1, P2, P3, P0]).transpose()
print('X (Spherical LTCS):\n{}'.format(X))

dH = -P1[1]

d1S = P1[0]
d2S = P2[0]
d3S = P3[0]
alpha = P2[1] - P1[1]
beta = P3[1] - P2[1]

P1p = numpy.array([data[0,0] * math.cos(data[0,1]),\
                  data[0,0] * math.sin(data[0,1]),\
                  data[0,2]])
P2p = numpy.array([data[1,0] * math.cos(data[1,1]),\
                  data[1,0] * math.sin(data[1,1]),\
                  data[1,2]])
P3p = numpy.array([data[2,0] * math.cos(data[2,1]),\
                  data[2,0] * math.sin(data[2,1]),\
                  data[2,2]])

delta12 = P2p - P1p
delta13 = P3p - P1p

d13 = math.sqrt(delta13[2]**2 + delta13[0]**2)
sigma13 = math.atan2(delta13[2], delta13[0])
sigma31 = math.atan2(delta13[0], delta13[2])
omega = math.pi - (alpha + beta)
sigma1c = sigma13 - beta
sigma3c = sigma31 + alpha
d1c = d13 * math.sin(alpha)/math.sin(omega)
d3c = d13 * math.sin(beta)/math.sin(omega)
zc = P1p[2] + d1c * math.sin(sigma1c)
print('zc: {}'.format(zc))
zc = P3p[2] + d3c * math.sin(sigma3c)
print('zc: {}'.format(zc))
xc = P1p[0] + d1c * math.cos(sigma1c)

Pc = numpy.array([xc, 0., zc])
deltac1 = P1p - Pc
sigmac1 = math.atan2(deltac1[2], deltac1[0])
deltac2 = P2p - Pc
sigmac2 = math.atan2(deltac2[2], deltac2[0])
psi = sigmac2
eta = 2*math.pi - (sigmac2 - sigmac1)
d1s = d13 * math.sin(eta)/math.sin(alpha + beta)
d3s = d13 * math.sin(psi)/math.sin(alpha + beta)
sigma1s = sigma13 + psi
sigma3s = sigma31 - eta
# xs = P1p[0] + d1s * math.cos(sigma1s)
xs = P3p[0] + d3s * math.cos(sigma3s)
zs = P1p[2] + d1s * math.cos(sigma3s)
ys = P1p[1] + dH
Ps = numpy.array([xs, ys, zs])
print('Ps: {}'.format(Ps))

Y = numpy.vstack([P1p, P2p, P3p, Ps]).transpose()
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