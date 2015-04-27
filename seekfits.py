#!/usr/bin/python
# M. Bogosavljevic, AOB, March 2015
# expanding capability of images.py by Z. Ioannou 
# *** work in progress ****
# ---------------------------------------------------------------
# This listens to a directory for presence of new fits files
# reads the header of the files, calculates the mean value for fun and prints. 
# It will work until KeyboardInterrupt, which is CTRL-C
#
# Run similar to this 
# python seekfits.py /mypath/ timeintervalinseconds 
# python seekfits.py ./test/ 3
#

import os, time
import pyfits
import numpy as np

def seekfits(path_to_watch, tsleep):

  before = dict ([(f, None) for f in os.listdir (path_to_watch)])

  # counter needed just to print out progress
  count = 0.

  try:
    while 1:
         count = count + 1
         print int(count), "*LISTENING*", path_to_watch, time.strftime('%X %x %Z')
         after = dict ([(f, None) for f in os.listdir (path_to_watch)])
         added = [f for f in after if not f in before]

         if added: 
           print "Added: ", ", ".join (added)
           for filein in enumerate(added):                  
             # check if it is a fits file
             filename = path_to_watch + filein[1]          
             # WARNING - hardcoded the '.fits' or '.fit' extensions
             if (filename.endswith('.fits') or filename.endswith('.fit')):
               # Can load both data and header with this trick
               data, hdr = pyfits.getdata(filename, header=True)              

               # Get Observation Header for image and display it
               dateobs     = hdr['DATE-OBS']
               timeobs     = hdr['UTC']
               exposure    = hdr['EXPTIME']
               CCDfilter   = hdr['FILTER']

         before = after
         time.sleep (tsleep)   # Wait for tsleep seconds before repeating
  except KeyboardInterrupt:
    pass

if __name__ == "__main__":

    import sys

    # check if path to watch ends with '/'
    path_to_watch = sys.argv[1]
    if not(path_to_watch.endswith('/')):
      path_to_watch = path_to_watch + '/'

    # this is the time period in seconds between checkin the path directory
    tsleep  = float(sys.argv[2])

    seekfits(path_to_watch, tsleep)
