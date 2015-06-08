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
!!!! UNFINISHED !!!!
"""
import os, time
import glob
import PIL
#import Image
import matplotlib.pyplot as plt

#===============================================================================
# Select last N frames of FITS files based on a filename wildcard.
#===============================================================================
def makefilelist(path,wildcard,nframes):

    # Move to the specified directory
    os.chdir(path)
    # Search for all files in the current directory that satisfy the wildcard provided
    filelist = sorted(glob.glob(wildcard))
    # print ("I see: ",filelist)
    lastn    = filelist[-nframes:]
    return lastn

#===============================================================================
# Main program
#===============================================================================
def loop_gif(ref_filename,wildcard,nframes,tsleep):

    # convert to int if string
    nframes = int(nframes)
    tsleep  = int(tsleep)

    path, filename = os.path.split(ref_filename)    
    png_dir = path + '/png'

    lastn    = makefilelist(png_dir,wildcard,nframes)
    print ('Lastn',lastn)
    oneimage = plt.imread(lastn[0]) 
    mylist   = []
    mylist.append(oneimage)
    plt.show(block=False)

    try:
        while 1:
             for filein in lastn[-nframes-1:]:                  
                 im = plt.imread(filein)
                 mylist.append(im)
                 frame = plt.imshow(im)
                 plt.show(block=False)
                 time.sleep(1)

             time.sleep(tsleep)   # Wait for tsleep seconds before repeating

    except KeyboardInterrupt:
        pass

if __name__ == "__main__":

   import sys

   ref_filename    = sys.argv[1]
   wildcard        = sys.argv[2]
   nframes         = sys.argv[3]
   tsleep          = sys.argv[4]

   loop_gif(ref_filename,wildcard,nframes,tsleep)
