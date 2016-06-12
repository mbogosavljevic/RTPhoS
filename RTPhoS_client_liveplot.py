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

def minutes_before_now(sometime):
# check how many minutes sincce time of data taken (assuming UTC)
#  expects sometime as 1990-01-01|00:00:00.00
    now = datetime.utcnow()
    sometime2 = datetime.strptime(sometime, "%Y-%m-%d|%H:%M:%S.%f")
    elapsedTime = sometime2 - now
    minutes = elapsedTime.total_seconds() / 60.
    return minutes


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
    if args.port is not None:
      print "!RTPhoS: you have specified both port and logfile to monitor!"
      print "This is not possible. Continuing with using the logfile."
      uselog = True
else:
    uselog = False
    if args.port is None: 
        print "Set -port or -logfile!"
        break

# bandpass filter must be at the beginning of the messages sent
if args.filter is not None:
    print "RTPhoS: filtering messages, bandpass %s" % args.filter
    bandpass_filter = "{\"bandpass\": \"" + args.filter
else:
    print "-filter option not set!"
    break

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

#################################
# if using live updating log file
################################

if uselog:

######
     try: 
        # open the log file
        # read all data that is currently there
        with open(args.logile) as f:
          obsid         = zip(*[line.split() for line in f])[0]
          serverport    = zip(*[line.split() for line in f])[1]
          sendtimeUTC   = zip(*[line.split() for line in f])[2]
          bandpass      = zip(*[line.split() for line in f])[3]
          UTCdatetime   = zip(*[line.split() for line in f])[4]
          BJD           = zip(*[line.split() for line in f])[5] 
          targetflux    = zip(*[line.split() for line in f])[6]
          targetfluxerr = zip(*[line.split() for line in f])[7]
          compflux      = zip(*[line.split() for line in f])[8]
          compfluxerr   = zip(*[line.split() for line in f])[9]
          seeing        = zip(*[line.split() for line in f])[10]
        
        # remember the time when file has been read last
        # and length of file
        nline = len(obsid)
        lastdatafiletime = os.stat(args.logfile).st_mtime

        # TBD:
        # use just data that matches -filter argument
        # something like matches = [x for x in a if x=='a']

        # TBD:
        # now filter out the last -npoints
        # .....

        # convert UTCdatetime in minutes before now
        # make also other lists in the same loop
        for k in range(0,len(UTCdatetime)+1):
            newx = minutes_before_now(UTCdatetime[k])
            xdata = xdata.append(newx)
            yrawtarget=yrawtarget.append(targetflux[k])
            yrawtargeterr=yrawtargeterr.append(targetfluxerr[k])
            yrawcomp=yrawcomp.append(compflux[k])
            yrawcomperr=yrawcomperr.append(compfluxerr[k])
            # for now use the same raw counts for diff flux:
            ydflux=yrawtarget
            ydfluxerr=yrawtargeterr
            yseeing=yseeing.append(seeing[k])
            
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

        # Differential Photometry
        ax4.grid(True)
        ax4.margins(0.1,0.1)
        ax4.set_ylabel('Relative Flux')
        line4, = ax4.plot(xdata, yrawtarget, 'bo')
        fig.legend([line4,], ['Diff'], loc='upper left', \
                   bbox_to_anchor=[0.9,0.48], shadow=True, numpoints=1, prop={'size':8})

        # Seeing 
        ax5.grid(True)
        ax5.margins(0.1,0.1)
        ax5.set_xlabel('Time JD+ ') #+str(stripjd))
        ax5.set_ylabel('Pixels')
        line5, = ax5.plot(xdata, yseeing, 'bo') 
        fig.legend([line5,], ['Seeing'], loc='upper left', \
                   bbox_to_anchor=[0.9,0.27], shadow=True, numpoints=1, prop={'size':8})
######### while
        while True:
               # check if input data file has been modified
               datafiletime = os.stat(args.logfile).st_mtime
################
               if datafiletime != lastdatafiletime:
                   # update time of last file update
                   lastdatafiletime = datafiletime
                   # check the new length of the file 
                   filelen = file_len(args.targetfile)
                   print "File has now: ", filelen

####################
                   if nline <= filelen:
                      for row in range(nline,filelen+1):

                         # move to the line stated by row variable
                         with open(args.logfile,'r') as fp:
                            j = 1
                            for line in fp:
                             if j == row:
                              useline = line
                              break
                             else:
                              j = j + 1
                         
                         newobsid         = line[0]
                         newserverport    = line[1]
                         newsendtimeUTC   = line[2]
                         newbandpass      = line[3]
                         newUTCdatetime   = line[4]
                         newBJD           = line[5] 
                         newtargetflux    = line[6]
                         newtargetfluxerr = line[7]
                         newcompflux      = line[8]
                         newcompfluxerr   = line[9]
                         newseeing        = line[10]
                         newx = minutes_before_now(line[4])
 
                         xdata = xdata.append(newx)
                         yrawtarget=yrawtarget.append(newyrawtarget)
                         yrawtargeterr=yrawtargeterr.append(newyrawtargeterr)
                         yrawcomp=yrawcomp.append(newyrawcomp)
                         yrawcomperr=yrawcomperr.append(newyrawcomperr)
                         # for now use the same raw counts for diff flux:
                         ydflux=yrawtarget
                         ydfluxerr=yrawtargeterr
                         yseeing=yseeing.append(seeing[k])
                         
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
                         
                         ## Target thumbnail
                         # ax1.plot([14,14],[0,14],'r:')             # Plot cross-hairs
                         # ax1.plot([0,14],[14,14],'r:')             #      -""-
                         # ax1.imshow(target_crop, cmap='gray', norm=LogNorm(vmin=cropmin, vmax=maxintens))
                         
                         ## Comparison thumbnail
                         # ax2.plot([14,14],[0,14],'r:')             # Plot cross-hairs
                         # ax2.plot([0,14],[14,14],'r:')             #      -""-
                         # ax2.imshow(comp_crop, cmap='gray', norm=LogNorm(vmin=cropmin, vmax=maxintens))
                         
                         # Target Raw counts
                         line3a.set_xdata(xdata)
                         line3a.set_ydata(yrawtarget)
                         line3b.set_xdata(xdata)
                         line3b.set_ydata(yrawcomp)
                         ax3.plot([newx,newx],[newtargetflux-newtargetfluxerr, newtargetflux+newtargetfluxerr], 'r-')
                         ax3.plot([newx,newx],[newcompflux-newcompfluxerr, newcompflux+newcompfluxerr], 'g-')       
                         ax3.relim()
                         
                         # Differential flux
                         line4.set_xdata(xdata)
                         line4.set_ydata(yrawtarget)
                         ax4.plot([newx,newx],[newtargetflux-newtargetfluxerr, newtargetflux+newtargetfluxerr], 'b-')
                         ax4.relim()
                         
                         # Seeing plot
                         line5.set_xdata(xdata)
                         line5.set_ydata(yseeing)
                         ax5.relim()
                         
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
        break


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
