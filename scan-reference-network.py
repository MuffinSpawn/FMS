# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 16:05:49 2018

@author: peter
"""
import logging
import numpy
import numpy.linalg as linalg
import os
import signal
import sys
import platform
import time

from CESAPI.connection import *
from CESAPI.command import *
from CESAPI.packet import *

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def signal_handler(signal, frame):
        print('You pressed Ctrl+C!')
        sys.exit(0)

def press_any_key_to_continue():
    print('<<< Press any key to continue... >>>')
    if platform.system() == 'Windows':
        os.system('pause')
    else:  # assuming Linux
        os.system('read -s -n 1')

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
        logger.debug('Initializing...')
        command.Initialize()
    
    logger.debug('setting measurement mode...')
    command.SetMeasurementMode(ES_MM_Stationary)  # ES_MeasMode (only choice for AT4xx)

    logger.debug('setting stationary mode parameters...')
    mode_params = StationaryModeDataT()
    mode_params.lMeasTime = 1000  # 1 second
    command.SetStationaryModeParams(mode_params)

    logger.debug('setting coordinate system type...')
    command.SetCoordinateSystemType(ES_CS_SCC)  # one of ES_CoordinateSystemType

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
    settings.bSendReflectorPositionData = int(1)
    settings.bTryMeasurementMode = int(0)
    settings.bHasNivel = int(1)
    settings.bHasVideoCamera = int(1)
    command.SetSystemSettings(settings)

# TODO: implement Ciddor & Hill with IOF update trick
def ciddor_and_hill(env_params):
    return 0.0

def set_refraction_index(command):
    min_refraction_index = 1.000150
    max_refraction_index = 1.000331
    mid_refraction_index = (min_refraction_index + max_refraction_index)/2.0
    refraction_params = command.GetRefractionParams()
    if refraction_params.dIfmRefractionIndex <= mid_refraction_index:
        refraction_params.dIfmRefractionIndex = max_refraction_index
    else:
        refraction_params.dIfmRefractionIndex = min_refraction_index
    command.SetRefractionParams(refraction_params)

    env_params = command.GetEnvironmentParams()
    index_of_refraction = ciddor_and_hill(env_params)
    command.SetRefractionParams(refraction_params)

def measure(command, setiof=False):
        if setiof:
            set_refraction_index(command)

        return command.StartMeasurement()

def main():
    signal.signal(signal.SIGINT, signal_handler)
    connection = LTConnection()
    try:
        logger.info('Connecting to the laser tracker...')
        connection.connect()
        command = CommandSync(connection)

        print("Acquire an initialization reflector.")
        press_any_key_to_continue()
        logger.info('Initializing laser tracker...')
        initialize(command, manualiof=False)
        
        measurements = []
        
        ordinal_strings = ['first', 'second', 'third']
        for index in range(3):
            print("Acquire the {} reference reflector.".format(ordinal_strings[index]))
            press_any_key_to_continue()
            logger.info('Measuring reflector..')
            measurements.append(measure(command, setiof=False))

        init_coordinates_LTCS = numpy.ndarray((3, 3))
        for index,measurement in enumerate(measurements):
            init_coordinates_LTCS[index,0] = measurement.dVal1
            init_coordinates_LTCS[index,1] = measurement.dVal2
            init_coordinates_LTCS[index,2] = measurement.dVal3
        print('Initial LTCS Coordinates:\n{}'.format(init_coordinates_LTCS))
    finally:
        connection.disconnect()

if __name__ == '__main__':
    main()
