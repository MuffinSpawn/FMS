# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 13:38:45 2018

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
        stream = connection.connect()

        while (True):
          unread_count = stream.unreadCount()
          if unread_count > 0:
            in_packet = stream.read()
            packet_type = packetType(in_packet)
            logger.info('Packet Type: {}'.format(packet_type))
            time.sleep(0.2)


    finally:
        connection.disconnect()

if __name__ == '__main__':
    main()
