#!/usr/bin/python
# 
"""
M. Bogosavljevic, AOB, June 2015
This code greps the RTPhoS log output 
and plots a calibrated source flux

Inputs:
- log : filename for rtphos log output (updated in real time)
- source : string containing desired source name
- comp1  : string containing start to use as comparison [required]
- comp2  : [optional] must be "none" if skipped
- comp3  : [optional] must be "none" if skipped
- nlast  : start plotting with nlast measurements

Outputs:
- just the matplotlib plot
https://docs.python.org/2/library/argparse.html
"""
import os, time

def phot_monitor(log, source, comp1, comp2, comp3, nlast, tsleep):
    
    print "=== RTPhoS Photometry Monitor Started === " + time.strftime('%X %x %Z')
   
    fsource = []
    fcomp1  = []
    fcomp2  = []
    fcomp3  = []

    try:
        while 1:

            # this command I used when there were no object names in optphot
            # now it can be done more easily

            # grep FRAME_TIME
            # grep source
            # grep comp1, comp2, comp3 etc.

            mycmd = " grep -B 1 -A 7 \"FRAME_TIME\" " + log + \
                " | grep -v \"\\-\\-\" |  awk '{printf \"%s%s\", $0, (NR%8?OFS:RS)}' |" + \
                " sort -k2 | awk '{print $2 \" \" $5 \" \" $10 \" \" $11 \" \" $25 \" \" $26}' "
            
            print mycmd
            os.system(mycmd)

# TBD
            # WORK IN PROGRESS!

            time.sleep(tsleep)   # Wait for tsleep seconds before repeating
    except KeyboardInterrupt:
        pass

if __name__ == "__main__":

   import sys
   log             = sys.argv[1]
   source          = sys.argv[2]
   comp1           = sys.argv[3]
   comp2           = sys.argv[4]
   comp3           = sys.argv[5]
   nlast           = sys.argv[6]   
   tsleep          = int(sys.argv[7])

   phot_monitor(log, source, comp1, comp2, comp3, nlast, tsleep)

