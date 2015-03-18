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
             filename=filein[1]
             # WARNING - hardcoded the '.fits' extension
             if filename.endswith('.fits'):
               # Can load both data and header with this trick
               data, hdr=pyfits.getdata(filename, header=True)              
               # Get Observation Header for image and display it
               obsdate = hdr['DATE-OBS']
               print obsdate
               # Display the data mean
               print "Mean value:",np.mean(data)
           
         before = after
         time.sleep (tsleep)   # Wait for tsleep seconds before repeating

  except KeyboardInterrupt:
     pass

if __name__ == "__main__":

    import sys
    path_to_watch = sys.argv[1]
    tsleep  = float(sys.argv[2])
    seekfits(path_to_watch, tsleep)
