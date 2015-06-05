#!/usr/bin/python
# 
"""
M. Bogosavljevic, AOB, June 2015
this should listen to a directory (obtained from ref_filename)
for files containing namestring wildcard
and find the last nframes there
Check if there is a gif of the fits file
and then create an animated gif from those last nframes
**** Work in progress ****
"""

import astropy.io.fits as pyfits
import os
import datetime

def loop_gif(ref_filename,namestring,nframes):

   path, filename = os.path.split(ref_filename)    

   if not os.path.exists(path+'/gifs'):
       print "Creating directory"
       print path+'/gifs'

   before = dict ([(f, None) for f in os.listdir(path)])

   try:
      while 1:
          print "*GIF LOOP RUNNING* ", path, time.strftime('%X %x %Z')
          after = dict ([(f, None) for f in os.listdir(path)])
          added = [f for f in after if not f in before]

          # Get all the fits files into a file list.
          fitsfiles = makefilelist(namestring)
  
          if added: 
                  count = count + 1             
                  print "Added files: ", ", ".join (added)
                  for filein in added:                  
                     # check if it is a fits file
                     filename = dirs['data']+'/'+filein
                     print filename
                     # WARNING - hardcoded the '.fits' or '.fit' extensions
                     if (filename.endswith('.fits') or filename.endswith('.fit')):
                        # Can load both data and header with this trick
                        data2, hdr = pyfits.getdata(filename, header=True)              
                        print("I read image: "+filename)

                before = after
                time.sleep(tsleep)   # Wait for tsleep seconds before repeating

        except KeyboardInterrupt:
           pass

if __name__ == "__main__":

   import sys

   path    = sys.argv[1]
   nframes = sys.argv[2]

   loop_gif(ref_filename,path,nframes)
