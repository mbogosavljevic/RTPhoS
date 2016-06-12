#!/usr/bin/env python
# RTPhoS client module which receives data from server module
# and logs into a flat ascii file.

# Usage:
# RTPhoS_client_logger.py [port] 
# {plus optional params and switches, do -h or --help to get options}

# port should be a string like tcp://localhost:port or other address
# Examples:
# RTPhoS_client_logger.py tcp://localhost:5556 
# RTPhoS_client_logger.py tcp://92.96.59.99:5556 

import os
import sys
import argparse
import zmq
import time
from datetime import datetime
import numpy as np
from astropy.io import fits
import json
