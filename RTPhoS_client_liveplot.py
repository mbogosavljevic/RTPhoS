#!/usr/bin/env python
# RTPhoS live plotting client module which either
# a) receives data from server module live, in which case nothing is saved 
# or
# b) monitors and plots the updating of logged data created by RTPhoS_client_logger.py

# Usage:
# RTPhoS_client_liveplot.py {-port [port] | -logfile [filename]} -filter bandpass
# plus optional params and switches, do -h or --help to get options

# port should be a string like tcp://localhost:port or other address
# Examples:
# RTPhoS_client_liveplot.py -port tcp://localhost:5556 
# RTPhoS_client_liveplot.py -logfile mylogfile.dat 
# if both port and logfile are specified, the code will complain
# and continue with using the logfile, NOT port

import os
import sys
import argparse
import time
from datetime import datetime
import numpy as np
from astropy.io import fits
import zmq
import json
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import rcParams

####################
# minutes_before_now
def minutes_before_now(sometime):
# check how many minutes sincce time of data taken (assuming UTC)
#  expects sometime as 1990-01-01|00:00:00.00
    now = datetime.utcnow()
    sometime2 = datetime.strptime(sometime, "%Y-%m-%d|%H:%M:%S")
    elapsedTime = sometime2 - now
    minutes = elapsedTime.total_seconds() / 60.
    return minutes

#################
# file_len
# Returns the actual number of lines in a text file
# do not have lines with just whitespace!
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


##########################################################
### MAIN CODE ###
##########################################################

# parse the arguments from the command line
parser = argparse.ArgumentParser(description='RTPhoS client live plotter')

# mandatory params:
# 
parser.add_argument('-filter', type=str, \
                    help='If set, will only plot data mathing the bandpass given as string')
###
parser.add_argument('-port', metavar='port', type=str, \
                    help='TCP/IP port to listen')
# or
parser.add_argument('-logfile', type=str, \
                    help='Optional string, file name for the log to be saved')
####

# optional
parser.add_argument('-npoints', type=int, \
                    help='If reading logfile, also plot last npoints')

parser.add_argument('-bjd', type=str, \
                    help ='Use BJD as time axis')

# TBD:
#parser.add_argument('-pastmin', type=int, \
#                    help='If reading lofile, also plot last pastmin minutes of data')

args = parser.parse_args()

################

if args.logfile is not None:
    uselog = True
    if args.port is not None:
      print "!RTPhoS: you have specified both port and logfile to monitor!"
      print "This is not possible. Continuing with using the logfile."
else:
    uselog = False
    if args.port is None: 
        print "Set -port or -logfile!"


# bandpass filter must be at the beginning of the messages sent
if args.filter is not None:
    print "RTPhoS: filtering messages, bandpass %s" % args.filter
    bandpass_filter = "{\"bandpass\": \"" + args.filter
else:
    print "-filter option not set!"


bandpass_filter = "{\"bandpass\": \"" + args.filter

print "RTPhoS: Live plot with bandpass filter: ", args.filter
       
# Initialize data lists
xdata=[]         # X-axis data (time, either past minutes or BJD)
yrawtarget = []  # Raw target counts
yrawtargeterr=[] # Raw target error bars
yrawcomp=[]      # Raw comparison counts
yrawcomperr=[]   # Raw comparison error bars
ydflux=[]        # Differential photometry counts
ydfluxerr=[]     # Differential photometry error bars
yseeing=[]       # Seeing data
UTCtimes = []
#################################
# if using live updating log file
################################

if uselog:

######
     try: 
        # open the log file
        # read all data that is currently there
        nline = 0
        # nline has value 0 for line 1 in file
        with open(args.logfile,'r') as fp:
          for line in fp:
              line = line.strip()  
              columns = line.split()
              nline = nline + 1
              # skip comment line
              if columns[0] == '#RTPhoS:':
                  print "skipped:", columns
              else:
                  obsid         = columns[0]
                  serverport    = columns[1]
                  UTCdatetime   = columns[2]
                  bandpass      = columns[3]
                  BJD           = np.float64(columns[4]) 
                  targetflux    = np.float32(columns[5])
                  targetfluxerr = np.float32(columns[6])
                  compflux      = np.float32(columns[7])
                  compfluxerr   = np.float32(columns[8])
                  seeing        = np.float16(columns[9])

                  # TBD:
                  # use just data that matches -filter argument
                  # something like matches = [x for x in a if x=='a']
                  # UTCtimes is a list, stores all actual times for
                  # recaluclating minutes_before_now later
                  UTCtimes.append(UTCdatetime)
                  newx = minutes_before_now(UTCdatetime)
                  print UTCdatetime, round(newx,3)
                  xdata.append(round(newx,3))
                  yrawtarget.append(targetflux)
                  yrawtargeterr.append(targetfluxerr)
                  yrawcomp.append(compflux)
                  yrawcomperr.append(compfluxerr)
                  # for now use the same raw counts for diff flux:
                  ydflux=yrawtarget
                  ydfluxerr=yrawtargeterr
                  yseeing.append(seeing)

        # remember the time when file has been read last
        lastdatafiletime = os.stat(args.logfile).st_mtime

        # TBD:
        # now filter out the last -npoints
        # .....

        # Initialize Plotting Environment
        plt.ion()
        fig = plt.figure(figsize=(12,12))
        ax1  = fig.add_subplot(421)
        ax2  = fig.add_subplot(422)
        ax3  = fig.add_subplot(412)
        ax4  = fig.add_subplot(413, sharex=ax3)
        ax5  = fig.add_subplot(414, sharex=ax3)

        # make plot settings once
        ax1.text(1, 5, 'Target', color='yellow', fontsize=10)
        ax1.set_xlabel("Physical X: ")
        ax1.set_ylabel("Physical Y: ") 
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
         
        ax2.text(1, 5, 'Comp', color='yellow', fontsize=10)
        ax2.set_xlabel("Physical X: ") 
        ax2.set_ylabel("Physical Y: ")
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])

        # Target & Comparison raw counts
        ax3.grid(True)
        ax3.margins(0.1,0.1)
        ax3.set_ylabel('Raw Counts')
        line3a, = ax3.plot(xdata, yrawtarget, 'ro')
        line3b, = ax3.plot(xdata, yrawcomp, 'go') 
        fig.legend([line3a, line3b,], ['Target', 'Comp'], loc='upper left', \
                bbox_to_anchor=[0.9,0.68], shadow=True, numpoints=1, prop={'size':8})
        # plot the error bars
        errorlines3a = []
        errorlines3b = []
        for k in range(0,len(xdata)):
            l, = ax3.plot([xdata[k],xdata[k]],[yrawtarget[k]-yrawtargeterr[k],\
                                          yrawtarget[k]+yrawtargeterr[k]], 'r-')
            errorlines3a.append(l)
            l, = ax3.plot([xdata[k],xdata[k]],[yrawcomp[k]-yrawcomperr[k],\
                                          yrawcomp[k]+yrawcomperr[k]], 'r-')
            errorlines3b.append(l)

        # Differential Photometry
        ax4.grid(True)
        ax4.margins(0.1,0.1)
        ax4.set_ylabel('Relative Flux')
        line4, = ax4.plot(xdata, ydflux, 'bo')
        fig.legend([line4,], ['Diff'], loc='upper left', \
                   bbox_to_anchor=[0.9,0.48], shadow=True, numpoints=1, prop={'size':8})
        # plot the error bars and collect them in
        # a new variable
        errorlines4 = []
        for k in range(0,len(xdata)):
            l, = ax4.plot([xdata[k],xdata[k]],[ydflux[k]-ydfluxerr[k],\
                                          ydflux[k]+ydfluxerr[k]], 'r-')
            errorlines4.append(l)

        # Seeing 
        ax5.grid(True)
        ax5.margins(0.1,0.1)
        if args.bjd:
            ax5.set_xlabel('Time JD+ ') 
        else:
            ax5.set_xlabel('Minutes [now = 0]')
 
        ax5.set_ylabel('Pixels')
        line5, = ax5.plot(xdata, yseeing, 'bo') 
        fig.legend([line5,], ['Seeing'], loc='upper left', \
                   bbox_to_anchor=[0.9,0.27], shadow=True, numpoints=1, prop={'size':8})

 
        npoints_current = len(xdata)
        fig.canvas.draw()
######### while
        while True:
               # check if input data file has been modified
               datafiletime = os.stat(args.logfile).st_mtime
################
               
               if datafiletime != lastdatafiletime:
                   # update time of last file update
                   lastdatafiletime = datafiletime
                   # check the new length of the file 
                   filelen = file_len(args.logfile)
                   print "File has now: ", filelen
                   # recalculate xdata for old points

                   xdata = []
                   for t in UTCtimes:
                       newtime = minutes_before_now(t)
                       xdata.append(round(newtime,3)) 
                   
####################
                   print "Starting to read from line:", nline+1
                   if (nline+1) <= filelen:
                     for row in range(nline,filelen):
                         # move to the line stated by row variable
                         with open(args.logfile,'r') as fp:
                            j = 0
                            for line in fp:
                             if j == row:
                              useline = line
                              break
                             else:
                              j = j + 1
                         
                         useline = useline.strip()  
                         useline = useline.split()
                         newobsid         = useline[0]
                         newserverport    = useline[1]
                         newUTCdatetime   = useline[2]
                         newbandpass      = useline[3]
                         
                         # TBD:
                         # here we shoud make a check if the data is
                         # still in the same bandpass!
                         newBJD           = np.float64(useline[4]) 
                         newtargetflux    = np.float32(useline[5])
                         newtargetfluxerr = np.float32(useline[6])
                         newcompflux      = np.float32(useline[7])
                         newcompfluxerr   = np.float32(useline[8])
                         newseeing        = np.float16(useline[9])

                         UTCtimes.append(newUTCdatetime)

                         newx = minutes_before_now(newUTCdatetime)
                         xdata = np.append(xdata,round(newx,3))
                         yrawtarget=np.append(yrawtarget,newtargetflux)
                         yrawtargeterr=np.append(yrawtargeterr,newtargetfluxerr)
                         yrawcomp=np.append(yrawcomp,newcompflux)
                         yrawcomperr=np.append(yrawcomperr,newcompfluxerr)
                         # for now use the same raw counts for diff flux:
                         ydflux=yrawtarget
                         ydfluxerr=yrawtargeterr
                         yseeing=np.append(yseeing,newseeing)
                         
                         # ## TBD:
                         # medianintens = np.median(target_crop)
                         # target_crop[target_crop==0] = medianintens # Remove zero values
                         # target_crop = np.log(target_crop)          # Use for log scale plotting
                         # comp_crop[comp_crop==0] = medianintens     # Remove zero values
                         # comp_crop = np.log(comp_crop)              # Use for log scale plotting
                         # cropmin = np.amin(target_crop)
                         # cropmax = np.amax(target_crop)
                         ## Attempt for a reasonable intensity scale 
                         # maxintens = ((cropmax-cropmin)/2.0)+cropmin

                         # clear the plot
                         
                         ## Target thumbnail
                         # ax1.plot([14,14],[0,14],'r:')             # Plot cross-hairs
                         # ax1.plot([0,14],[14,14],'r:')             #      -""-
                         # ax1.imshow(target_crop, cmap='gray', norm=LogNorm(vmin=cropmin, vmax=maxintens))
                         
                         ## Comparison thumbnail
                         # ax2.plot([14,14],[0,14],'r:')             # Plot cross-hairs
                         # ax2.plot([0,14],[14,14],'r:')             #      -""-
                         # ax2.imshow(comp_crop, cmap='gray', norm=LogNorm(vmin=cropmin, vmax=maxintens))
                         
                         # Target Raw counts
                         #print len(xdata)
                         #rint len(yrawtarget)
                         #line3a, = ax3.plot(xdata, yrawtarget, 'ro')
                         #line3b, = ax3.plot(xdata, yrawcomp, 'go') 
                         
                         line3a.set_xdata(xdata)
                         line3a.set_ydata(yrawtarget)
                         line3b.set_xdata(xdata)
                         line3b.set_ydata(yrawcomp)
                         ax3.relim()
                         
                         safe_lines3a = errorlines3a
                         safe_lines3b = errorlines3b
                         errorlines3a = []
                         errorlines3b = []
                         # remove previous error bars
                         for k in range(0,npoints_current):
                             ax3.lines.remove(safe_lines3a[k])
                             ax3.lines.remove(safe_lines3b[k])
                         for k in range(0,len(xdata)):
                             l, = ax3.plot([xdata[k],xdata[k]], \
                                           [yrawtarget[k]-yrawtargeterr[k],\
                                            yrawtarget[k]+yrawtargeterr[k]], 'r-')
                             errorlines3a.append(l)
                             l, = ax3.plot([xdata[k],xdata[k]], \
                                           [yrawcomp[k]-yrawcomperr[k],\
                                            yrawcomp[k]+yrawcomperr[k]], 'r-')
                             errorlines3b.append(l)

                         # Differential flux
                         line4.set_xdata(xdata)
                         line4.set_ydata(yrawtarget)
                         ax4.relim()
                         # plot the error bars
                         safe_lines4 = errorlines4
                         errorlines4 = []
                         for k in range(0,npoints_current):
                             ax4.lines.remove(safe_lines4[k])
                         for k in range(0,len(xdata)):
                             l, = ax4.plot([xdata[k],xdata[k]],\
                                           [ydflux[k]-ydfluxerr[k],\
                                            ydflux[k]+ydfluxerr[k]], 'r-')
                             errorlines4.append(l)

                         # Seeing plot
                         line5.set_xdata(xdata)
                         line5.set_ydata(yseeing)
                         ax5.relim()

                         npoints_current = len(xdata)
                         fig.canvas.draw()
######################### end for rows
#################### end if nline <= filelen
                     nline = row + 1                  
################ end if datafiletime != 
               else: 
                   # Sleep 1 second before attempting to parse input data file(s) again
                   time.sleep(1)
######### end while
###### end try
     except (KeyboardInterrupt, SystemExit):
        print "Process aborted by user."

###############################
## if live plotting from server
###############################
### TBD:
#else:
#     # Socket to talk to server
#     context = zmq.Context()
#     socket = context.socket(zmq.SUB)
#     socket.connect ("%s" % port)
#     socket.setsockopt(zmq.SUBSCRIBE,bandpass_filter)
#     print "RTPhoS: Collecting updates from broadcasting server ", port
#
#     try: 
#        while True:
#              # get the message packet, which is a dictionary
#              messagedata =json.loads(socket.recv())
#
#              now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#        
#              obsid         = str(messagedata['obsid'])
#              serverport    = str(messagedata['port'])
#              sendtimeUTC   = str(messagedata['sendtimeUTC'])
#              bandpass      = str(messagedata['bandpass'])
#              UTCdatetime   = str(messagedata['UTCdatetime'])
#              BJD           = str(messagedata['BJD'])
#              targetflux    = str(messagedata['targetflux'])
#              targetfluxerr = str(messagedata['targetfluxerr'])
#              compflux      = str(messagedata['compflux'])
#              compfluxerr   = str(messagedata['compfluxerr'])
#              seeing        = str(messagedata['seeing'])
#      
#              timage       = np.array(messagedata['thumbnail1'])
#              cimage       = np.array(messagedata['thumbnail2'])
#
#              print "RTPhoS: message received %s" % now
#              print "Sender: %s" % obsid
#              print "Sent time UTC : %s" % sendtimeUTC
#              print "Observation time UTC: %s" % UTCdatetime
#              now = datetime.utcnow()
#              sometime = datetime.strptime(UTCdatetime, "%Y-%m-%d|%H:%M:%S.%f")
#              elapsedTime = sometime - now
#              secondold = elapsedTime.total_seconds()
#              print "Data is %s seconds old." % str(round(abs(secondold),1))
#              print "--- Message logged."
#     except (KeyboardInterrupt, SystemExit):
#        print "Process aborted by user."
#        break
