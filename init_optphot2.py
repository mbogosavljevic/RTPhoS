#!/usr/bin/python
####################################
## M. Bogosavljevic, SQU, March 2015
##### for use in the run_rtphos.py pipeline
##### WORK IN PROGRESS 
## INPUT
## Needs only the XPApoint name from DS9
## (you can also run it from shell, example ./init_optphot.py 7f000101:59109)
##
## OUTPUT
## psf.cat   - input list for Optphot PSF procedure
## stars.cat - input list for Optphot measurements 
## (by design, they differ only in the target source
##  which is the first source in stars.cat)
##
## Usage:
## 1 - start DS9 from your current working directory
## 2 - load ds9_optphot.ans in DS9 Analysis pull down menu
## 3 - load an image and create regions for several stars
## 4 - centroid regions in DS9 (Region -> Select All, Region -> Ceontrid)
## 5 - click on a region of the target so that it is the only one selected
## 6 - clik on "Opthot Initialize" in DS9 Analysis pull down menu
##

import pyregion
import ds9
import sys

def readregions(xpapoint):

    # create a ds9 object linked with an XPA point
    win = ds9.ds9(xpapoint)

    ## redo centering just in case - this did not work great - why?
    #  win.set("regions centroid radius 20 | iteration 10")
    #  print(" Regions centroided ")

    # initally, target should be the only selected region
    targetrows = win.get("regions selected") 
    target = pyregion.parse(targetrows)
    xt = target[0].coord_list[0]
    yt = target[0].coord_list[1]

    # check if target is length of one
    if len(target)!=1:
       print(" Selected one region (only)! ")
       raise Exception(" Selected one region (only)! ")
    else:
       print('%-15s %7.2f %7.2f' % ("Target (X,Y):", xt, yt) )
       
    # invert selection of regions (for comparison/psf stars)
    win.set("regions select invert")
    compstars = win.get("regions selected")
    # go back to the first selected region
    win.set("regions select invert")

    # comparison stars 
    comps = pyregion.parse(compstars)
   
    print('%-25s' % ("Comparison stars at:") )

    #pdb.set_trace()
    # Write out output to files that Optphot expects
    # these are psf.cat for objects to be usd for psf
    # and stars.cat for all objects to be measured
    text_file1 = open("psf.cat", "w")
    text_file2 = open("stars.cat", "w")
    text_file1.write("!\n!\n!\n")
    text_file2.write("!\n!\n!\n")
 
    # write the target as the first source in stars.cat
    text_file2.write('%-5s %8.1f %8.1f \n' % ("1", xt, yt) )
 
    ncomp = len(comps)
    for k in range(0,ncomp):
      x = comps[k].coord_list[0]
      y = comps[k].coord_list[1]
      text_file1.write('%-5i %8.1f %8.1f \n' % (k+1, x, y) )
      text_file2.write('%-5i %8.1f %8.1f \n' % (k+2, x, y) )
      print('%-5i %8.1f %8.1f' % (k+1, x, y) )
 
    text_file1.close()
    text_file2.close()

if __name__ == "__main__":

    import sys
    xpaname = sys.argv[1]
    print("Using xpaname: ", xpaname)
    readregions(xpaname)
