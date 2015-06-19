#!/usr/bin/python
# M. Bogosavljevic, AOB, June 2015
# based on wuclient.py from http://zguide.zeromq.org/py:all
#
#   Photometry update client
#   Connects SUB socket to tcp://[inputIP]:5556

import time

def minutes_passed(oldepoch, mins):
    return time.time() - oldepoch >= mins * 60

def phot_client(soruceIP, objectID):

    import zmq

    #  Socket to talk to server
    context = zmq.Context()
    socket = context.socket(zmq.SUB)
    
    # zmq requires convert to unicode string
    objectID = objectID.decode('ascii')

    print("Collecting updates from phot server - start - "+time.strftime('%X %x %Z'))
 
#   socket.connect("tcp://"+sourceIP+":5556")
    # if localhost, otherwise uncomment above
    socket.connect("tcp://localhost:5556")
    socket.setsockopt_string(zmq.SUBSCRIBE, objectID)

    starttime = time.time()
    # Listen until CTRL-C
    try:
        while 1:
          message = socket.recv_string()
          uttime, filt, mag, magerr = string.split()
          print (uttime, filt, mag, magerr, time.strftime('%X %x %Z'))
          # remind that you are listening every 1 minute
          if minutes_passed(starttime,1):
              print (" Still listening "+ time.strftime('%X %x %Z'))
    except KeyboardInterrupt:
        pass
       

if  __name__ == "__main__":

    import sys
    sourceIP       = sys.argv[1]
    objectID       = sys.argv[2]
    phot_client(sourceIP, objectID)

