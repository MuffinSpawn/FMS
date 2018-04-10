# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 16:05:49 2018

@author: Peter G. Lane (petergwinlane@gmail.com)
"""

import csv
import logging
import math
import numpy
import numpy.linalg as linalg
import os
import signal
import sys
import platform
import time

import CESAPI.connection
import CESAPI.command
from CESAPI.packet import *
import CESAPI.refract

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Convert Spherical Counter-Clockwise coordinates to Right-Handed Rectangular coordinates
def convert_network_to_RHR(coordinates_SCC):
    cordinates_RHR = numpy.ndarray(numpy.shape(coordinates_SCC))
    for index,point in enumerate(coordinates_SCC):
        cordinates_RHR[index,0] = point[0]*math.cos(point[1])*math.sin(point[2])
        cordinates_RHR[index,1] = point[0]*math.sin(point[1])*math.sin(point[2])
        cordinates_RHR[index,2] = point[0]*math.cos(point[2])
    return cordinates_RHR

# Apply the transform matrix to Cartesian coordinates
def transform_network(transform_matrix, coordinates_RHR_in):
    X = numpy.vstack([numpy.ones(numpy.shape(coordinates_RHR_in)[0]), coordinates_RHR_in.transpose()])
    return numpy.dot(transform_matrix, X)[1:,:].transpose()

# Convert Right-Handed Rectangular coordinates to Cylindrical Counter-Clockwise coordinates
def convert_network_to_CCC(coordinates_RHR):
    cordinates_CCC = numpy.ndarray(numpy.shape(coordinates_RHR))
    for index,cartesian_point in enumerate(coordinates_RHR):
        cordinates_CCC[index,0] = math.sqrt(numpy.sum(cartesian_point[:2]**2))
        cordinates_CCC[index,1] = math.atan2(cartesian_point[1], cartesian_point[0])
        cordinates_CCC[index,2] = cartesian_point[2]
    return cordinates_CCC

def main(coordinates_string):
    transform_matrix_string = '''
0.9999999999999998,-2.0513105103425744e-16,-4.209414934674971e-16,-7.796286064106449e-16
1347.5240143985754,-0.949485899932389,-0.31380098496329817,0.002338304197711084
408.1015988490014,0.002529933335232659,-0.00020343901273492393,0.9999967790194383
1986.5095112337217,-0.3137994985140675,0.9494887574096931,0.0009870580448882466
    '''
    transform_matrix = numpy.array(list(map(lambda x: x.split(','), transform_matrix_string.split('\n')[1:-1])), dtype=numpy.float32)

    coordinates_SCC = numpy.array([coordinates_string.split(','),], dtype=numpy.float32)
    coordinates_CCC = convert_network_to_CCC(transform_network(transform_matrix, convert_network_to_RHR(coordinates_SCC)))

    print('LTCS Coordinates: {}'.format(coordinates_SCC))
    print('DSCS Coordinates: {}'.format(coordinates_CCC))
    
if __name__ == '__main__':
    coordinates_string = '2206.847993,  -0.871582,  1.318256'
    coordinates_string = '2539.36, -0.666054, 1.29189'
    coordinates_string = '1924.73, -1.375285, 1.839'   # Table
    coordinates_string = '2058.4,0.928388,2.21644'     # Room corner by door
    main(coordinates_string)
