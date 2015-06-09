#!/usr/bin/python
# 
"""
M. Bogosavljevic, AOB, June 2015
this should listen to a directory (obtained from ref_filename)
subdirecrtory /png 
for images containing namestring wildcard
and find the last nframes there
and then create an animated gif from those last nframes
"""
import os, time
import glob

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
    mylist   = ' '.join(lastn)
    outfile = 'last'+str(nframes)+'loop.gif'  
    outmade = os.path.isfile(outfile)

    try:
        while 1:
            newlastn    = makefilelist(png_dir,wildcard,nframes)
            newlist   = ' '.join(newlastn)
            print newlist
            if (newlist != mylist) or (not(outmade)):
                mycmd = 'convert -delay 50 -loop 0 '+newlist+' '+outfile
                print mycmd
                os.system(mycmd)
                mylist = newlist
                outmade = True
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
