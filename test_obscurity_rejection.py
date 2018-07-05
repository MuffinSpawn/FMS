# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 14:27:45 2018

@author: plane
"""

import math
import os
import os.path as path
import platform
import sys
import numpy
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.collections as mcollections
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt  # Matplotlib's pyplot: MATLAB-like syntax

angles = numpy.arange(-math.pi, math.pi, 0.1)
radii = numpy.array([500, 1000, 1500, 2000, 2500, 3000])

fig = plt.figure()
plt.plot(0, 0, 'b+')
for radius in radii:
    xs = radius*numpy.cos(angles)
    ys = radius*numpy.sin(angles)
    plt.plot(xs, ys, 'r.')
plt.axes().set_aspect('equal', 'datalim')

theta0 = -math.pi/2.0+1.0
#theta0 = -math.pi
#theta0 = math.pi
print(theta0)
r0 = 1000
width = 800
depth = 1000
m0 = math.tan(theta0 + math.pi/2.0)
c = math.tan(math.pi/2.0 - theta0)/math.cos(theta0)
b1 = (r0 - depth/2.0) * c
b2 = (r0 + depth/2.0) * c
for radius in radii:
    print(math.acos(width/(2*radius)))
    theta_r_over_2 = math.pi/2.0 - math.acos(width/(2*radius))
    print('Theta0: {}\tTheta_r/2: {}'.format(theta0, theta_r_over_2))
    theta1 = theta0 - theta_r_over_2
    plt.plot([0,radius*math.cos(theta1)], [0,radius*math.sin(theta1)], 'b-')
    theta2 = theta0 + theta_r_over_2
    print('Theta1: {}\tTheta2: {}'.format(theta1, theta2))
    plt.plot([0,radius*math.cos(theta2)], [0,radius*math.sin(theta2)], 'b-')
    no_plot_condition = None
    if theta1 < -math.pi:
        theta1 += 2*math.pi
        no_plot_condition = numpy.logical_or(angles>theta1, angles<theta2)
    elif theta2 > math.pi:
        theta2 -= 2*math.pi
        no_plot_condition = numpy.logical_or(angles>theta1, angles<theta2)
    else:
        no_plot_condition = numpy.logical_and(angles>theta1, angles<theta2)
    #plot_angles = angles[numpy.logical_not(numpy.logical_and(numpy.logical_and(no_plot_condition, radius >= r1), radius <= r2))]
    no_plot_angles = angles[no_plot_condition]
    print('No Plot Angles: {}'.format(no_plot_angles))
    ms = numpy.tan(no_plot_angles)
    print('Slopes: {}'.format(ms))
    x1s = b1 / (ms - m0)
    y1s = ms * x1s
    r1s = numpy.sqrt(x1s**2 + y1s**2)
    print('R1s: {}'.format(r1s))
    x2s = b2 / (ms - m0)
    y2s = ms * x2s
    r2s = numpy.sqrt(x2s**2 + y2s**2)
    print('R2s: {}'.format(r2s))
    yes_plot_angles = no_plot_angles[numpy.logical_or(radius < r1s, radius > r2s)]
    plot_angles = numpy.hstack((angles[numpy.logical_not(no_plot_condition)], yes_plot_angles))
    xs = radius*numpy.cos(plot_angles)
    ys = radius*numpy.sin(plot_angles)
    plt.plot(xs, ys, 'g.')
