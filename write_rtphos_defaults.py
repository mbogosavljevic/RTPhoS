#!/usr/bin/python
# M. Bogosavljevic, AOB, June 2015

def write_rtphos_defaults( pathdefs, pathdata, pathbias, pathdark, pathflat, stringbias, stringdark, stringflat, namebias, namedark, nameflat, filejpl, fileleap, cprefix, sradius, aradius, cradius, starnumber, skyskew, skyfit, gain, verbose):
 
    defs_file = open(pathdefs, "w")

    defs_file.write(pathdefs   + "     # File containing new settings\n")
    defs_file.write(pathdata   + "     # Raw data frames\n")		
    defs_file.write(pathbias   + "     # Bias frames for calibration\n")	
    defs_file.write(pathdark   + "     # Dark frames for calibration\n")
    defs_file.write(pathflat   + "     # Flat frames for calibration\n")	
    defs_file.write(stringbias + "     # bias frames wildcard\n")		
    defs_file.write(stringdark + "     # dark frames wildcard\n")		
    defs_file.write(stringflat + "     # flat frames wildcard\n")		
    defs_file.write(namebias   + "     # Masterbias filename\n")		
    defs_file.write(namedark   + "     # Masterdark filename\n")		
    defs_file.write(nameflat   + "     # Masterflat filename\n")		
    defs_file.write(filejpl    + "     # path to JPL Ephemeris file\n")	     
    defs_file.write(fileleap   + "     # path to leap seconds data file\n")
    defs_file.write(cprefix    + "     # Prefix for calibrated files\n")
    defs_file.write(sradius    + "     # Search radius for optimal centroiding\n")	
    defs_file.write(aradius    + "     # Aperature radius for photometry\n")	
    defs_file.write(cradius    + "     # PSF Clipping radius for optimal photom\n")	
    defs_file.write(starnumber + "     # Star number to optimize photometry for\n")	
    defs_file.write(skyskew    + "     # Sky profile skew for optimal photometry\n")	
    defs_file.write(skyfit     + "     # Sky fitting switch for optimal photometry\n")	
    defs_file.write(gain       + "     # Instrument Gain e-/ADU\n")	
    defs_file.write(verbose    + "     # Verbose output switch\n")	

    defs_file.close()

if  __name__ == "__main__":

    import sys
    pathdefs       = sys.argv[1]
    pathdata       = sys.argv[2]
    pathbias       = sys.argv[3]
    pathdark       = sys.argv[4]
    pathflat       = sys.argv[5]
    stringbias     = sys.argv[6]
    stringdark     = sys.argv[7]
    stringflat     = sys.argv[8]
    namebias       = sys.argv[9]
    namedark       = sys.argv[10]
    nameflat       = sys.argv[11]
    filejpl        = sys.argv[12]
    fileleap       = sys.argv[13]
    cprefix        = sys.argv[14]
    sradius        = sys.argv[15]
    aradius        = sys.argv[16]
    cradius        = sys.argv[17]
    starnumber     = sys.argv[18]
    skyskew        = sys.argv[19]
    skyfit         = sys.argv[20]
    gain           = sys.argv[21]
    verbose        = sys.argv[22]

    write_rtphos_defaults( pathdefs, pathdata, pathbias, pathdark, pathflat, stringbias, stringdark, stringflat, namebias, namedark, nameflat, filejpl, fileleap, cprefix, sradius, aradius, cradius, starnumber, skyskew, skyfit, gain, verbose)
