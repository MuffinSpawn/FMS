# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 11:52:46 2018

@author: plane
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

def main():
    signal.signal(signal.SIGINT, signal_handler)

    connection = CESAPI.connection.Connection()
    try:
        logger.info('Connecting to the laser tracker...')
        connection.connect()
        command = CESAPI.command.CommandSync(connection)

        logger.info('Getting compensations...')
        compensation = command.GetCompensations2()
        print('{} compensations available.'.format(compensation.iTotalCompensations))
        '''
        self.iTotalCompensations = int(0)
        self.iInternalCompensationId = int(0)
        self.cTrackerCompensationName = str()  # 32 bytes max
        self.cTrackerCompensationComment = str()  # 128 bytes max
        self.cADMCompensationName = str()  # 32 bytes max
        self.cADMCompensationComment = str()  # 128 bytes max
        self.bHasMeasurementCameraMounted = int(0)
        self.bIsActive = int(0)
        '''
    finally:
        connection.disconnect()

if __name__ == '__main__':
    main()
