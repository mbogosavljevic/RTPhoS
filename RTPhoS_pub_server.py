import zmq
import sys
import time
import random
from astropy.io import fits
import json
import numpy as np
from collections import OrderedDict

# get command line parameters port, obsid, bandpass, thumbnail1
# the rest will be generated random numbers
# example command line:
# python pub_server.py 5556 Vidoje R test.fits

# for current development these are command line arguments
port =  sys.argv[1]
obsid = sys.argv[2]
bandpass = sys.argv[3]
hdulist1 = fits.open(sys.argv[4])
scidata1 = hdulist1[0].data

# connect to a publishing port
context = zmq.Context()
socket = context.socket(zmq.PUB)
socket.bind("tcp://*:%s" % port)

localUTtime = 0.

try:
    while True:
        # convert numpy array for stamp to list
        scidata_list1 = scidata1.tolist()

        # create some random numbers for time, flux etc
        localUTtime = localUTtime + 1.0
        BJD = random.randrange(0, 101)
        targetflux = random.randrange(0,100000)
        targetfluxerr = np.sqrt(targetflux)
        compflux = random.randrange(0,200000)
        compfluxerr = np.sqrt(compflux)
        seeing = float("{0:.2f}".format(random.uniform(1,2)))

        # for fun print to terminal
        print "-----------------------------"
        print "bandpass", bandpass
        print "obsid", obsid
        print "port", port
        print "localUTtime", localUTtime
        print "BJD", BJD
        print "targetflux", targetflux
        print "targetfluxerr", targetfluxerr
        print "compflux", compflux
        print "compfluxerr", compfluxerr
        print "seeing", seeing

        # create an ordered dictionary
        message = OrderedDict ( [("bandpass", bandpass), ("obsid", obsid),  ("port", port), \
                                 ("localUTtime", localUTtime), ("BJD",BJD), \
                                 ("targetflux", targetflux), ("targetfluxerr", targetfluxerr),\
                                 ("compflux", compflux), ("compfluxerr", compfluxerr),\
                                 ("seeing", seeing),  ("thumbnail1",scidata_list1) ] )

        # convert to json message 
        jsonmessage = json.dumps(message)
        socket.send(jsonmessage)
        
        # debug
        #print(jsonmessage)

        # for now, send every 1 second, hardcoded
        time.sleep(1)
except (KeyboardInterrupt, SystemExit):
    print "Process aborted by user."
    pass
