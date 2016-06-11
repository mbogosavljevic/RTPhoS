#!/usr/bin/env python
# RTPhoS server module which watches output photometry files
# and broadcasts new data over TCP/IP using ZeroMQ

# Usage:
# RTPhoS_pub_server.py [port] [obsid] [band] [targetfile] [tsleep] 
# {plus optional params and switches}

import argparse
import zmq
import sys
import time
import os
from astropy.io import fits
import json
import numpy as np
from collections import OrderedDict
from datetime import datetime

#########################################################################################

def minutes_before_now(sometime):
# check how many minutes sincce time of data taken (assuming UTC)
#  expects sometime as 1990-01-01|00:00:00
    now = datetime.utcnow()
    sometime2 = datetime.strptime(sometime, "%Y-%m-%d|%H:%M:%S")
    elapsedTime = sometime2 - now
    minutes = elapsedTime.total_seconds() / 60.
    return minutes

#########################################################################################

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

#########################################################################################

def datatext_to_message(textfile, lastlineread, firsttime):
# Used to read data or comparison star data files with the format as:
#   No.    Date                BJD           time_err  Flux     Flux_err    seeing flag  filename
#    1  1990-01-01|00:00:00  2447892.500058   5.00   33368.6875 647.047546   5.46 O gauss01.fits 
    with open(textfile,'r') as fp:
         i = -1
         if firsttime:
             for line in fp:
                 i = i + 1
                 save_i = i
             useline = line
             print useline
         else:
             for line in fp:
                 i = i + 1
                 if i == lastlineread:
                     useline = line
                     save_i = i

         last = useline
         line = last.strip()   
         print "Here:", line
         columns = line.split()
         thisentry         = columns[0]
         UTCdatetime  = columns[1]
         diffminutes  = minutes_before_now(UTCdatetime)
         BJD       = columns[2] 
         flux      = columns[4]
         fluxerr   = columns[5]
         seeing    = columns[6]

         # create a dictionary
         part_message = {"lastread":save_i, "dUTCminutes":diffminutes, "BJD":BJD, \
                                "flux":flux, "fluxerr":fluxerr, "seeing":seeing} 
         return part_message

#########################################################################################

def do_the_work(args,t_part_message,row,firsttime):
    if args.thumbtarget is not None:
        print "Reading target fits file: %s" % args.thumbtarget
        hdulist1 = fits.open(args.thumbtarget)
        timage = hdulist1[0].data
        timage_list = timage.tolist()
    else:
        timage_list = np.array(float('NaN')).tolist()
        if args.compfile is not None:
            c_part_message = datatext_to_message(args.compfile, row, firsttime)
        if args.thumbcomp is not None:
            print "Reading comp fits file: %s" % args.thumbcomp
            hdulist2 = fits.open(args.thumbcomp)
            cimage = hdulist2[0].data
        else:
            cimage_list = np.array(float('NaN')).tolist()

    dUTCminutes =  t_part_message['dUTCminutes']
    BJD = t_part_message['BJD']
    targetflux = t_part_message['flux']
    targetfluxerr = t_part_message['fluxerr']
    seeing = t_part_message['seeing']

    print "-----------------------------"
    print "Read data entry:", t_part_message['lastread']+1
    print "bandpass", args.band
    print "obsid", args.obsid
    print "port", args.port
    print "dUTCminutes", dUTCminutes
    print "BJD", BJD
    print "targetflux", targetflux
    print "targetfluxerr", targetfluxerr
    print "seeing", seeing
    if compdone:
        compflux = c_part_message['flux']
        compfluxerr = c_part_message['fluxerr']
        print "compflux", compflux
        print "compfluxerr", compfluxerr
    else:
        compflux = float('NaN')
        compfluxerr = float('NaN')

    
    # create an ordered dictionary
    message = OrderedDict ( [("bandpass", args.band), ("obsid", args.obsid),  ("port", args.port), \
                  ("dUTCminutes", dUTCminutes), ("BJD",BJD), ("targetflux", targetflux), \
                  ("targetfluxerr", targetfluxerr), ("compflux", compflux), \
                  ("compfluxerr", compfluxerr), ("seeing", seeing), \
                  ("thumbnail1",timage_list), ("thumbnail2",cimage_list)] )
 
    # convert to json message 
    jsonmessage = json.dumps(message)
    socket.send(jsonmessage)
    print jsonmessage
   
##########################################################
# MAIN ###################################################
#########################################################

# parse the arguments from the command line
parser = argparse.ArgumentParser(description='RTPhoS live data broadcasting module')

# required params
parser.add_argument('port', metavar='port', type=int, \
                    help='TCP/IP port to be used for the broadcast')
parser.add_argument('obsid', metavar='obsid', type=str, \
                    help='Observatory ID string to be used')
parser.add_argument('band', metavar='band', type=str, \
                    help='String stating which photometric band is being broadcast')
parser.add_argument('targetfile', metavar='targetfile', type=str, \
                    help='Target object data file to be monitored for updates and broadcast.')
parser.add_argument('tsleep', metavar='tsleep', type=float, \
                    help='Sleep time between checks of target data file')

# boolean switch
parser.add_argument('-public', action='store_true', \
                    help='Set this switch to broadcast the target flux. Default FALSE')
# optional
parser.add_argument('--compfile', metavar='compfile', type=str, \
                    help='Comparison star data to be monitored for updates and broadcast.')
parser.add_argument('--thumbtarget', metavar='thumbtarget', type=str, \
                    help='50x50 pix fits thumbnail of the target object')
parser.add_argument('--thumbcomp', metavar='thumbcomp', type=str, \
                    help='50x50 pix fits thumbnail of the comparison star')
parser.add_argument('--path', metavar='path', type=str, \
                    help='path to data files to be read [default: current dir]')

parser.set_defaults(path='./')
args = parser.parse_args()

# connect to a publishing port
context = zmq.Context()
socket = context.socket(zmq.PUB)
socket.bind("tcp://*:%s" % args.port)

lastdatafiletime = 0
lastdataentry = -1
firsttime = True 

try:
  while True:
      # check if input data file has been modified
      datafiletime = os.stat(args.targetfile).st_mtime
      compdone = False
      if datafiletime != lastdatafiletime:
          # update time when update noted
          lastdatafiletime = datafiletime
          # check the new length of the file 
          filelen = file_len(args.targetfile)
          print "I see the file already has", filelen, " lines"
          # read the just last line if it is the first time
          # otherwise, read all lines that are new
          if firsttime:
              print lastdataentry
              t_part_message = datatext_to_message(args.targetfile, lastdataentry, firsttime)
              do_the_work(args,t_part_message,lastdataentry,firsttime)
              lastdataentry = t_part_message['lastread']
              firsttime = False
          else:
              for row in range(lastdataentry,filelen):
                  t_part_message = datatext_to_message(args.targetfile, row, firsttime)
                  do_the_work(args,t_part_message,row,firsttime)
              lastdataentry = t_part_message['lastread']+2
      else: 
          # Sleep before attempting to parse input data file(s) again
          time.sleep(args.tsleep)
          print "Checking again...", datetime.now()

except (KeyboardInterrupt, SystemExit):
  print "Process aborted by user."
  pass
