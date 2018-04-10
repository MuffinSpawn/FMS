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

def signal_handler(signal, frame):
        print('You pressed Ctrl+C!')
        sys.exit(0)

def press_any_key_to_continue():
    if platform.system() == 'Windows':
        os.system('pause')
    else:  # assuming Linux
        print('<<< Press any key to continue... >>>')
        os.system('read -s -n 1')
    print()

def initialize(command, forceinit=False, manualiof=False):
    units = SystemUnitsDataT()
    units.lenUnitType = ES_LU_Millimeter  # ES_LengthUnit
    # units.angUnitType = ES_AU_Radian  # ES_AngleUnit
    # units.tempUnitType = ES_TU_Celsius  # ES_TemperatureUnit
    # units.pressUnitType = ES_PU_Mbar  # ES_PressureUnit
    # units.humUnitType = ES_HU_RH  # ES_HumidityUnit
    logger.debug('Setting units...')
    command.SetUnits(units)
    
    status = command.GetSystemStatus()
    logger.debug('Tracker Processor Status: {}'.format(status.trackerProcessorStatus))
    if forceinit or status.trackerProcessorStatus != ES_TPS_Initialized:  # ES_TrackerProcessorStatus
        try:
            init_result = command.Initialize()
            logger.debug('Initialize status: {}'.format(init_result.packetInfo.status))
        except Exception as e:
            logger.error('Initialize failed: {}'.format(e))
            return
        # At least the AT401 seems to complain about an unknown command failing due to "the sensor" not being stable
        # on the next command after an initialize. The tracker is fine after that, so just ignore this as a bug in the firmware.
        try:
            status = command.GetSystemStatus()
            logger.debug('Tracker Processor Status: {}'.format(status.trackerProcessorStatus))
        except Exception as e:
            if not 'Command 64' in str(e):
                raise e
    
    logger.debug('setting measurement mode...')
    command.SetMeasurementMode(ES_MM_Stationary)  # ES_MeasMode (only choice for AT4xx)

    logger.debug('setting stationary mode parameters...')
    mode_params = StationaryModeDataT()
    mode_params.lMeasTime = 1000  # 1 second
    command.SetStationaryModeParams(mode_params)

    logger.debug('setting coordinate system type to Right-Handed Rectangular...')
    command.SetCoordinateSystemType(ES_CS_RHR)  # one of ES_CoordinateSystemType

    logger.debug('setting system settings...')
    settings = SystemSettingsDataT()
    # one of ES_WeatherMonitorStatus
    if manualiof:
        settings.weatherMonitorStatus = ES_WMS_ReadOnly
    else:
        settings.weatherMonitorStatus = ES_WMS_ReadAndCalculateRefractions
    settings.bApplyStationOrientationParams = int(1)
    settings.bKeepLastPosition = int(1)
    settings.bSendUnsolicitedMessages = int(1)
    settings.bSendReflectorPositionData = int(0)
    settings.bTryMeasurementMode = int(0)
    settings.bHasNivel = int(1)
    settings.bHasVideoCamera = int(1)
    command.SetSystemSettings(settings)

def measure(command, rialg=None):
        CESAPI.refract.SetRefractionIndex(command, rialg)
        return command.StartMeasurement()

def loadDSCS(filename):
    with open(filename, 'r') as file:
        string_data = list(csv.reader(file, delimiter=','))

    reflector_names = []
    cylindrical_reflector_coordinates = numpy.ndarray((len(string_data), 3))
    for index,row in enumerate(string_data):
        reflector_names.append(row[0])
        point = list(map(lambda x: float(x), row[1:]))
        cylindrical_reflector_coordinates[index,:] = point

    cartesian_reflector_coordinates = numpy.ndarray((len(string_data), 3))
    for index,point in enumerate(cylindrical_reflector_coordinates):
        cartesian_reflector_coordinates[index,0] = point[0] * math.cos(point[1])
        cartesian_reflector_coordinates[index,1] = point[0] * math.sin(point[1])
        cartesian_reflector_coordinates[index,2] = point[2]
    return (numpy.array(reflector_names), cylindrical_reflector_coordinates, cartesian_reflector_coordinates)

def measure_plane_reflectors(command, rialg=CESAPI.refract.RI_ALG_Leica):
    refraction_index_algorithm = CESAPI.refract.AlgorithmFactory().refractionIndexAlgorithm(CESAPI.refract.RI_ALG_Leica)
    measurements = []
    reflector_names = []

    command.FindReflector(1000)

    ordinal_strings = ['first', 'second', 'third']
    for index in range(3):
        reflector_name = input(\
            "\nEnter the name of the {} plane reflector:  ".format(ordinal_strings[index]))
        reflector_names.append(reflector_name)
        
        print("\nAcquire the {} reference reflector.".format(ordinal_strings[index]))
        press_any_key_to_continue()
        logger.info('Measuring reflector..')
        measurements.append(measure(command, rialg=refraction_index_algorithm))

    coordinates_LTCS = numpy.ndarray((3, 3))
    for index,measurement in enumerate(measurements):
        coordinates_LTCS[index,0] = measurement.dVal1
        coordinates_LTCS[index,1] = measurement.dVal2
        coordinates_LTCS[index,2] = measurement.dVal3
    return (reflector_names, coordinates_LTCS)

# Calculate the transform between three points from a source CS
# to three points from a target CS.
def calculate_transform(coordinates_RHR_in, coordinates_RHR_out):
    # calculate the tracker position (S') in the output CS to use as the fourth point
    logger.debug('Input Cartesian coordinates:\n{}'.format(coordinates_RHR_in))
    logger.debug('Output Cartesian coordinates:\n{}'.format(coordinates_RHR_out))
    A = coordinates_RHR_in[0]
    B = coordinates_RHR_in[1]
    C = coordinates_RHR_in[2]
    S = numpy.array([0, 0, 0])
    AB = B-A
    AC = C-A
    AS = -A
    z = numpy.cross(AB, AC)
    z_hat = z / linalg.norm(z)
    logger.debug('z_hat: {}'.format(z_hat))
    a = numpy.dot(z_hat, AS)
    b = numpy.dot(AB, AS)
    c = numpy.dot(AC, AS)

    Ap = coordinates_RHR_out[0]
    Bp = coordinates_RHR_out[1]
    Cp = coordinates_RHR_out[2]
    ABp = Bp-Ap
    ACp = Cp-Ap
    zp = numpy.cross(ABp, ACp)
    zp_hat = zp / linalg.norm(zp)
    logger.debug("z'_hat: {}".format(zp_hat))
    M = numpy.vstack([zp_hat, ABp, ACp])
    M_inv = linalg.inv(M)
    ASp = numpy.dot(M_inv, numpy.array([a, b, c]))
    logger.debug("AS': {}".format(ASp))
    Sp = ASp + Ap
    logger.debug("S': {}".format(Sp))
    
    # calculate the LTCS-to-DSCS transform matrix
    X = numpy.vstack([[1, 1, 1, 1], numpy.vstack([A, B, C, S]).transpose()])
    Y = numpy.vstack([[1, 1, 1, 1], numpy.vstack([Ap, Bp, Cp, Sp]).transpose()])
    return numpy.dot(Y, linalg.pinv(X))

def calculate_approx_LTCS(cartesian_DSCS, transform_matrix):
    X = numpy.vstack([[1, 1, 1, 1, 1], cartesian_DSCS.transpose()])
    return numpy.dot(transform_matrix, X)[1:,:].transpose()

def measurement_to_array(measurement):
    measurement_array = numpy.ndarray(9)
    measurement_array[0] = measurement.dVal1
    measurement_array[1] = measurement.dVal2
    measurement_array[2] = measurement.dVal3
    measurement_array[3] = measurement.dStd1
    measurement_array[4] = measurement.dStd2
    measurement_array[5] = measurement.dStd3
    measurement_array[6] = measurement.dTemperature
    measurement_array[7] = measurement.dPressure
    measurement_array[8] = measurement.dHumidity
    return measurement_array

def scan_reference_network(command, cartesian_LTCS):
    # adjust laser tracker
    logger.debug('setting coordinate system type to Counter-Clockwise Spherical...')
    command.SetCoordinateSystemType(ES_CS_SCC)  # one of ES_CoordinateSystemType
    
    logger.debug('Cartesian LTCS:\n{}'.format(cartesian_LTCS))
    spherical_points = numpy.ndarray(numpy.shape(cartesian_LTCS))
    for index,cartesian_point in enumerate(cartesian_LTCS):
        spherical_points[index,2] = math.sqrt(numpy.sum(cartesian_point**2))
        spherical_points[index,0] = math.atan2(cartesian_point[1], cartesian_point[0])
        spherical_points[index,1] = math.acos(cartesian_point[2]/spherical_points[index,2])
    logger.debug('Spherical LTCS:\n{}'.format(spherical_points))

    logger.info('Measuring reference network with (both faces)...')
    measurements_face1 = numpy.ndarray((numpy.shape(cartesian_LTCS)[0], 9))
    measurements_face2 = numpy.ndarray((numpy.shape(cartesian_LTCS)[0], 9))
    for index,spherical_point in enumerate(spherical_points):
        logger.debug('Directing laser to coordinates {}...'.format(spherical_point))
        command.GoPosition(int(1), spherical_point[0], spherical_point[1], spherical_point[2])

        # the tracker always switches to face 1 after a GoPosition command
        measurements_face1[index] = measurement_to_array(measure(command))

        # the tracker always switches to face 1 after a GoPosition command
        command.ChangeFace()
        measurements_face2[index] = measurement_to_array(measure(command))
    logger.debug('Face 2 Measurements:\n{}'.format(measurements_face2))
    logger.debug('Face 1 Measurements:\n{}'.format(measurements_face1))

    # calculate the average coordinates from the two-face measurements
    spherical_LTCS = numpy.ndarray((numpy.shape(cartesian_LTCS)))
    for index in range(numpy.shape(measurements_face1)[0]):
        spherical_LTCS[index] = (measurements_face1[index,:3] + measurements_face2[index,:3]) / 2.0
    return spherical_LTCS

# Apply the transform matrix to Cartesian coordinates
def transform_network(transform_matrix, coordinates_RHR_in):
    X = numpy.vstack([numpy.ones(numpy.shape(coordinates_RHR_in)[0]), coordinates_RHR_in.transpose()])
    return numpy.dot(transform_matrix, X)[1:,:].transpose()

# Convert Right-Handed Rectangular coordinates to Spherical Counter-Clockwise coordinates
def convert_network_to_SCC(coordinates_RHR):
    coordinates_SCC = numpy.ndarray(numpy.shape(coordinates_RHR))
    for index,cartesian_point in enumerate(coordinates_RHR):
        coordinates_SCC[index,2] = math.sqrt(numpy.sum(cartesian_point**2))  # r
        coordinates_SCC[index,0] = math.atan2(cartesian_point[1], cartesian_point[0])  # theta (azimuthal)
        coordinates_SCC[index,1] = math.acos(cartesian_point[2]/coordinates_SCC[index,2])  # phi (polar)
    return coordinates_SCC

# Convert Right-Handed Rectangular coordinates to Cylindrical Counter-Clockwise coordinates
def convert_network_to_CCC(coordinates_RHR):
    cordinates_CCC = numpy.ndarray(numpy.shape(coordinates_RHR))
    for index,cartesian_point in enumerate(coordinates_RHR):
        cordinates_CCC[index,0] = math.sqrt(numpy.sum(cartesian_point[:2]**2))
        cordinates_CCC[index,1] = math.atan2(cartesian_point[1], cartesian_point[0])
        cordinates_CCC[index,2] = cartesian_point[2]
    return cordinates_CCC

def print_configuration(transform_matrix, ref_spherical_LTCS, ref_cylindrical_DSCS, prop_spherical_LTCS, ds_spherical_LTCS, offset=1):
    for point_index,point in enumerate(ref_spherical_LTCS):
        for coordinate_index,coordinate in zip([1, 2, 0], point):
            if coordinate_index == 2:
                print('PredictedLTCSCoordinateSets[{},{}] = {:.3f}'.format(point_index+offset, coordinate_index, coordinate))
            else:
                print('PredictedLTCSCoordinateSets[{},{}] = {:.6f}'.format(point_index+offset, coordinate_index, coordinate))

    prop_offset = offset + numpy.shape(ref_spherical_LTCS)[0] + 1
    for point_index,point in enumerate(prop_spherical_LTCS):
        for coordinate_index,coordinate in zip([1, 2, 0], point):
            if coordinate_index == 2:
                print('PredictedLTCSCoordinateSets[{},{}] = {:.3f}'.format(point_index+prop_offset, coordinate_index, coordinate))
            else:
                print('PredictedLTCSCoordinateSets[{},{}] = {:.6f}'.format(point_index+prop_offset, coordinate_index, coordinate))

    ds_offset = offset + numpy.shape(ref_spherical_LTCS)[0] + 9
    for point_index,point in enumerate(ds_spherical_LTCS):
        for coordinate_index,coordinate in zip([1, 2, 0], point):
            if coordinate_index == 2:
                print('PredictedLTCSCoordinateSets[{},{}] = {:.3f}'.format(point_index+ds_offset, coordinate_index, coordinate))
            else:
                print('PredictedLTCSCoordinateSets[{},{}] = {:.6f}'.format(point_index+ds_offset, coordinate_index, coordinate))

    for point_index,point in enumerate(ref_cylindrical_DSCS):
        for coordinate_index,coordinate in zip([1, 2, 0], point):
            if coordinate_index == 2:
                print('MeasuredDSCSCoordinateSets[{},{}] = {:.3f}'.format(point_index+offset, coordinate_index, coordinate))
            else:
                print('MeasuredDSCSCoordinateSets[{},{}] = {:.6f}'.format(point_index+offset, coordinate_index, coordinate))

    print('TransformMatrix = NULL')
    transform_matrix_LTCS_DSCS = linalg.inv(transform_matrix)
    for row_index,row in enumerate(transform_matrix_LTCS_DSCS):
        for element_index,element in enumerate(row):
            print('TransformMatrix[{},{}] = {:.6f}'.format(row_index, element_index, element))

def define_unit_vectors(plane_coordinates_LTCS):
    P12 = plane_coordinates_LTCS[1] - plane_coordinates_LTCS[0]
    P13 = plane_coordinates_LTCS[2] - plane_coordinates_LTCS[0]
    z_hatp = P12/linalg.norm(P12)
    P13_x_z_hatp = numpy.cross(P13, z_hatp)
    x_hatp = P13_x_z_hatp/linalg.norm(P13_x_z_hatp)
    y_hatp = numpy.cross(z_hatp, x_hatp)
    return (x_hatp, y_hatp, z_hatp)

def calculate_ref_coordinates(plane_coordinates_LTCS, y_hatp, z_hatp):
    ref_coordinates = numpy.zeros(numpy.shape(plane_coordinates_LTCS))
    ref_coordinates[1,2] = numpy.dot(plane_coordinates_LTCS[1]-plane_coordinates_LTCS[0], z_hatp)
    ref_coordinates[2,1] = numpy.dot(plane_coordinates_LTCS[2]-plane_coordinates_LTCS[0], y_hatp)
    ref_coordinates[2,2] = numpy.dot(plane_coordinates_LTCS[2]-plane_coordinates_LTCS[0], z_hatp)
    return ref_coordinates

def calculate_Euler_angles(R):
    theta_x = math.atan2(R[2,1], R[2,2]) * 180 / math.pi
    theta_y = math.atan2(-R[2,0], math.sqrt(R[2,1]**2+R[2,2]**2)) * 180 / math.pi
    theta_z = math.atan2(R[1,0], R[0,0]) * 180 / math.pi
    return (theta_x, theta_y, theta_z)

def test():
    reflector_names = ['REF1', 'REF2', 'REF3']
    '''
    plane_coordinates_RHR = numpy.array([[ 1361.06082467,   172.89210929,  -411.51114447],
                                         [ 1901.78827007, -1463.23106064,  -413.21200562],
                                         [ 1903.76888187, -1460.93517541,   699.13857103]])
    '''
    plane_coordinates_RHR = numpy.array([[ 1901.78827007, -1463.23106064,  -413.21200562],
                                         [ 1361.06082467,   172.89210929,  -411.51114447],
                                         [ 1903.76888187, -1460.93517541,   699.13857103]])

    logger.debug('Y-Z Plane Cartesian coordinates:\n{}'.format(plane_coordinates_RHR))

    x_hatp, y_hatp, z_hatp = define_unit_vectors(plane_coordinates_RHR)
    logger.debug('Cartesian unit vectors:\n{}'.format(numpy.vstack([x_hatp, y_hatp, z_hatp])))

    # For confirmation (but mostly for fun), physically point out the three coordinate axes.
    '''
    connection = CESAPI.connection.Connection()
    try:
        logger.info('Connecting to the laser tracker...')
        connection.connect()
        command = CESAPI.command.CommandSync(connection)
        command.PointLaser(x_hatp[0], x_hatp[1], x_hatp[2])
        command.PointLaser(y_hatp[0], y_hatp[1], y_hatp[2])
        command.PointLaser(z_hatp[0], z_hatp[1], z_hatp[2])
        command.GoPosition(int(0), plane_coordinates_RHR[0,0], plane_coordinates_RHR[0,1], plane_coordinates_RHR[0,2])
    finally:
        connection.disconnect()

    other_point = numpy.ones((4,1))
    try:
        logger.info('Connecting to the laser tracker...')
        connection.connect()
        command = CESAPI.command.CommandSync(connection)
        logger.info('Initializing laser tracker...')
        initialize(command, manualiof=False)
        refraction_index_algorithm = CESAPI.refract.AlgorithmFactory().refractionIndexAlgorithm(CESAPI.refract.RI_ALG_Leica)
        logger.info('Measuring reflector..')
        measurement = measurement_to_array(measure(command, rialg=refraction_index_algorithm))[:3]
    finally:
        connection.disconnect()
    other_point[1,0] = measurement[0]
    other_point[2,0] = measurement[1]
    other_point[3,0] = measurement[2]
    '''

    ref_coordinates_RHR = calculate_ref_coordinates(plane_coordinates_RHR, y_hatp, z_hatp)
    logger.debug('Reference Cartesian coordinates:\n{}'.format(ref_coordinates_RHR))

    logger.info('Input Cartesian coordinates:\n{}'.format(plane_coordinates_RHR))
    logger.info('Output Cartesian coordinates:\n{}'.format(ref_coordinates_RHR))

    transform_matrix = calculate_transform(plane_coordinates_RHR, ref_coordinates_RHR)
    csv_output = ['LTCS-to-DSCS Transform Matrix:\n']
    for index,row in enumerate(transform_matrix):
        csv_output.append('{},{},{},{}\n'.format(row[0], row[1], row[2], row[3]))
    logger.info(''.join(csv_output))

    inv_transform_matrix = linalg.inv(transform_matrix)
    csv_output = ['DSCS-to-LTCS Transform Matrix:\n']
    for index,row in enumerate(inv_transform_matrix):
        csv_output.append('{},{},{},{}\n'.format(row[0], row[1], row[2], row[3]))
    logger.info(''.join(csv_output))

    tracker_coordinates_RHR = numpy.dot(transform_matrix, [[1.0],[0],[0],[0]])
    logger.info('Calculated tracker Cartesian coordinates:\n{}\n'.format(tracker_coordinates_RHR))    

    theta_x, theta_y, theta_z = calculate_Euler_angles(transform_matrix[1:,1:])
    logger.info('Euler angles:\n{}, {}, {}\n'.format(theta_x, theta_y, theta_z))

    ref_coordinates_SCC = convert_network_to_SCC(plane_coordinates_RHR)
    csv_output = ['Reference Spherical coordinates:\n']
    for index,point in enumerate(ref_coordinates_SCC):
        csv_output.append('{},{},{},{}\n'.format(reflector_names[index], point[2], point[0], point[1]))
    logger.info(''.join(csv_output))

    ref_coordinates_CCC = convert_network_to_CCC(ref_coordinates_RHR)
    csv_output = ['Reference Cylindrical coordinates (defined):\n']
    for index,point in enumerate(ref_coordinates_CCC):
        csv_output.append('{},{},{},{}\n'.format(reflector_names[index], point[0], point[1], point[2]))
    logger.info(''.join(csv_output))

    ref_coordinates_CCC = convert_network_to_CCC(numpy.dot(transform_matrix, numpy.vstack(([1, 1, 1], plane_coordinates_RHR.transpose())))[1:,:].transpose())
    csv_output = ['Reference Cylindrical coordinates (calculated):\n']
    for index,point in enumerate(ref_coordinates_CCC):
        csv_output.append('{},{:.3f},{:.6f},{:.3f}\n'.format(reflector_names[index], point[0], point[1], point[2]))
    logger.info(''.join(csv_output))

def main():
    signal.signal(signal.SIGINT, signal_handler)

    connection = CESAPI.connection.Connection()
    try:
        logger.info('Connecting to the laser tracker...')
        connection.connect()
        command = CESAPI.command.CommandSync(connection)
        
        print("Acquire an initialization reflector.")
        press_any_key_to_continue()
        logger.info('Initializing laser tracker...')
        initialize(command, forceinit=False, manualiof=False)

        reflector_names, plane_coordinates_RHR = measure_plane_reflectors(command)
        logger.debug('Y-Z Plane Cartesian coordinates:\n{}'.format(plane_coordinates_RHR.transpose()))

        x_hatp, y_hatp, z_hatp = define_unit_vectors(plane_coordinates_RHR)
        logger.debug('Cartesian unit vectors:\n{}'.format(numpy.vstack([x_hatp, y_hatp, z_hatp]).transpose()))

        # For confirmation (but mostly for fun), physically point out the three coordinate axes.
        command.PointLaser(x_hatp[0], x_hatp[1], x_hatp[2])
        command.PointLaser(y_hatp[0], y_hatp[1], y_hatp[2])
        command.PointLaser(z_hatp[0], z_hatp[1], z_hatp[2])
        command.GoPosition(int(0), plane_coordinates_RHR[0,0], plane_coordinates_RHR[0,1], plane_coordinates_RHR[0,2])
    finally:
        connection.disconnect()

    ref_coordinates_RHR = calculate_ref_coordinates(plane_coordinates_RHR, y_hatp, z_hatp)
    logger.debug('Reference Cartesian coordinates:\n{}'.format(ref_coordinates_RHR))

    logger.info('Input Cartesian coordinates:\n{}'.format(plane_coordinates_RHR))
    logger.info('Output Cartesian coordinates:\n{}'.format(ref_coordinates_RHR))

    transform_matrix = calculate_transform(plane_coordinates_RHR, ref_coordinates_RHR)

    csv_output = ['Transform Matrix:\n']
    for index,row in enumerate(transform_matrix):
        csv_output.append('{},{},{},{}\n'.format(row[0], row[1], row[2], row[3]))
    logger.info(''.join(csv_output))

    tracker_coordinates_RHR = numpy.dot(transform_matrix, [[1.0],[0],[0],[0]])
    logger.info('Calculated tracker Cartesian coordinates:\n{}\n'.format(tracker_coordinates_RHR))    

    theta_x, theta_y, theta_z = calculate_Euler_angles(transform_matrix[1:,1:])
    logger.info('Euler angles:\n{}, {}, {}\n'.format(theta_x, theta_y, theta_z))

    ref_coordinates_SCC = convert_network_to_SCC(plane_coordinates_RHR)
    csv_output = ['Reference Spherical coordinates:\n']
    for index,point in enumerate(ref_coordinates_SCC):
        csv_output.append('{},{},{},{}\n'.format(reflector_names[index], point[2], point[0], point[1]))
    logger.info(''.join(csv_output))

    ref_coordinates_CCC = convert_network_to_CCC(ref_coordinates_RHR)
    csv_output = ['Reference Cylindrical coordinates:\n']
    for index,point in enumerate(ref_coordinates_CCC):
        csv_output.append('{},{},{},{}\n'.format(reflector_names[index], point[0], point[1], point[2]))
    logger.info(''.join(csv_output))

if __name__ == '__main__':
    #main()
    test()
