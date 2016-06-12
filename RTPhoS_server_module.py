#!/usr/bin/env python
# RTPhoS server module which watches output photometry files
# and broadcasts new data as it arrives over TCP/IP using ZeroMQ

# Usage:
# RTPhoS_server_module.py [port] [obsid] [band] [targetfile] [tsleep] 
# {plus optional params and switches, do -h or --help to get options}

import os
import sys
import argparse
import zmq
import time
from datetime import datetime
import numpy as np
from astropy.io import fits
import numpy as np
import json
from collections import OrderedDict

#########################################################################################
# file_len
# Returns the actual number of lines in a text file
# do not have lines with just whitespace!
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

#########################################################################################
# datatext_to_message
# Parses the chosen line of the text file and returns a dictionary
# the convention for nline is that lines are numbered from 1 to file_len(textfile)
# Used to read target data or comparison star data files with the format as:
#   No.    UTCdatetime       BJD            terr[s]  Flux       Flux_err    seeing flag  filename
#    1  1990-01-01|00:00:00  2447892.500058   5.00   33368.6875 647.047546   5.46   O     gauss01.fits 
def datatext_to_message(textfile, nline):

    with open(textfile,'r') as fp:
         i = 1
         for line in fp:
             if i == nline:
                 useline = line
                 break
             else:
                 i = i + 1

         useline = useline.strip()   
         columns = useline.split()
         thisentry    = columns[0]
         UTCdatetime  = columns[1]
         BJD          = columns[2] 
         flux         = columns[4]
         fluxerr      = columns[5]
         seeing       = columns[6]

         # create a dictionary
         part_message = {"linenumber":i, "UTCdatetime":UTCdatetime, "BJD":BJD, \
                                "flux":flux, "fluxerr":fluxerr, "seeing":seeing} 
         return part_message

#########################################################################################
# format_and_broadcast
# Using the command line arguments in args, 
# check what parts of the message should be sent, read the given line of datafiles, 
# create json message and broadcast.
# For example, checks if thumbnails of the target and comparison exist and reads them in
# otherwise, it will broadcast NaN for the above
# checks if the comparison star data exists, etc.
# uses the datatext_to_message function defined above

def format_for_broadcast(args, t_part_message, nline):
# nline here serves only to find the same line in comparison star data file
# if that is the preference set 

    compdone = False
    if args.thumbtarget is not None:
        print "Reading target fits file: %s" % args.thumbtarget
        hdulist1 = fits.open(args.thumbtarget)
        timage = hdulist1[0].data
        timage_list = timage.tolist()
    else:
        timage_list = 'NaN'
        if args.compfile is not None:
            c_part_message = datatext_to_message(args.compfile, nline)
            compdone = True
        if args.thumbcomp is not None:
            print "Reading comp fits file: %s" % args.thumbcomp
            hdulist2 = fits.open(args.thumbcomp)
            cimage = hdulist2[0].data
        else:
            cimage_list = 'NaN'

    UTCdatetime =  t_part_message['UTCdatetime']
    BJD = t_part_message['BJD']
    if args.public:
        targetflux = t_part_message['flux']
        targetfluxerr = t_part_message['fluxerr']
    else:
        targetflux = 'NaN'
        # make target flux error a percentage in this case
        targetfluxerr = float(t_part_message['fluxerr']) / float(t_part_message['flux']) * 100.

    seeing = t_part_message['seeing']

    if compdone:
        compflux = c_part_message['flux']
        compfluxerr = c_part_message['fluxerr']
    else:
        compflux = 'NaN'
        compfluxerr = 'NaN'

    sendtimeUTC = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
    # create an ordered dictionary
    message = OrderedDict ( [("bandpass", args.band), ("sendtimeUTC", sendtimeUTC), ("obsid", args.obsid), \
                             ("port", args.port), \
                             ("UTCdatetime", UTCdatetime), ("BJD",BJD), ("targetflux", targetflux), \
                             ("targetfluxerr", targetfluxerr), ("compflux", compflux), \
                             ("compfluxerr", compfluxerr), ("seeing", seeing), \
                             ("thumbnail1", timage_list), ("thumbnail2", cimage_list)] )
 
    # convert to json message and broadcast 
    jsonmessage = json.dumps(message)
    return jsonmessage

##########################################################
### MAIN CODE ###
##########################################################

# parse the arguments from the command line
parser = argparse.ArgumentParser(description='RTPhoS live data broadcasting module')

# required params
parser.add_argument('port', metavar='port', type=str, \
                    help='TCP/IP port to be used for the broadcast')
parser.add_argument('obsid', metavar='obsid', type=str, \
                    help='Observatory ID string to be used')
parser.add_argument('band', metavar='band', type=str, \
                    help='String stating which photometric band is being broadcast')
parser.add_argument('targetfile', metavar='targetfile', type=str, \
                    help='Target object data file to be monitored for updates and broadcast.')
parser.add_argument('tsleep', metavar='tsleep', type=float, \
                    help='Sleep time between checks of target data file')

# boolean switch: set it in order to broadcast the fluxes, otherwise they will be private
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
port = int(args.port)
context = zmq.Context()
socket = context.socket(zmq.PUB)
socket.bind("tcp://*:%s" % port)

lastdatafiletime = 0
nline = 1
firsttime = True 

## main loop ##
try:
  while True:
      # check if input data file has been modified
      datafiletime = os.stat(args.targetfile).st_mtime
      if datafiletime != lastdatafiletime:
          # update time of last file update
          lastdatafiletime = datafiletime
          # check the new length of the file 
          filelen = file_len(args.targetfile)
          print "File has", filelen
          if firsttime:
              # read the just last line if it is the first time
              # otherwise, read all lines that are new since last visit
              print "RTPhoS: the input file has", filelen, " existing lines"
              nline = filelen
              t_part_message = datatext_to_message(args.targetfile, nline)
              jsonmessage = format_for_broadcast(args,t_part_message,nline)
              socket.send(jsonmessage)
              print "RTPhoS: Message sent ", datetime.now()
              print jsonmessage
              firsttime = False
              # move counter to next line for next read
              nline = nline + 1
          else:
              if nline <= filelen:
                  for row in range(nline,filelen+1):
                      t_part_message = datatext_to_message(args.targetfile, row)
                      jsonmessage = format_for_broadcast(args,t_part_message,row)
                      socket.send(jsonmessage)
                      now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                      print "RTPhoS: Message sent ", now
                      print jsonmessage
                      # set nline to one line more for next read
                      nline = row + 1                  
      else: 
          # Sleep before attempting to parse input data file(s) again
          time.sleep(args.tsleep)
          now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
          print "RTPhoS: Checking input files again...", now
except (KeyboardInterrupt, SystemExit):
  print "Process aborted by user."
  pass
