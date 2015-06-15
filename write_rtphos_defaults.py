#!/usr/bin/python
# M. Bogosavljevic, AOB, June 2015

def write_rtphos_defaults( rtphos, pathdefs, pathlog, pathans, pathdata, pathbias, pathdark, pathflat, stringbias, stringdark, stringflat, namebias, namedark, nameflat, cprefix, sradius, aradius, cradius, starnumber, skyskew, skyfit, gain, verbose, wildcard, nframes, tsleep):

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
    print "Wrote :"+pathdefs

    defs_file = open(pathans, "w")
    
    defs_file.write("param rtphosdef\n")
    defs_file.write(" rtphos  entry {RTPhoS location}       \""        +   rtphos + "\"\n" )
    defs_file.write(" pathdefs    entry {RTPhoS settings} \""              +   pathdefs   + "\"\n" )
    defs_file.write(" pathlog     entry {RTPhoS logfile}               \"" +   pathlog    + "\"\n" )
    defs_file.write(" pathans     entry {DS9 Analysis settings} \""        +   pathans    + "\"\n" )
    defs_file.write(" pathdata    entry {Raw data frames}	       \"" +   pathdata   + "\"\n" )   
    defs_file.write(" pathbias    entry {Bias frames for calibration}  \"" +   pathbias   + "\"\n" )    
    defs_file.write(" pathdark    entry {Dark frames for calibration}  \"" +   pathdark   + "\"\n" )    
    defs_file.write(" pathflat    entry {Flat frames for calibration}  \"" +   pathflat   + "\"\n" )   
    defs_file.write(" stringbias  entry {bias frames wildcard}	       \"" +   stringbias + "\"\n" )  
    defs_file.write(" stringdark  entry {dark frames wildcard}	       \"" +   stringdark + "\"\n" )  
    defs_file.write(" stringflat  entry {flat frames wildcard}	       \"" +   stringflat + "\"\n" )  
    defs_file.write(" namebias    entry {Masterbias filename}	       \"" +   namebias   + "\"\n" )  
    defs_file.write(" namedark    entry {Masterdark filename}	       \"" +   namedark   + "\"\n" )       
    defs_file.write(" nameflat    entry {Masterflat filename}	       \"" +   nameflat   + "\"\n" )
    defs_file.write(" cprefix     entry {Prefix for calibrated files}   \"" +  cprefix    + "\"\n" )
    defs_file.write(" sradius     entry {Search radius for optimal centroiding}     "   +  sradius + "\n" )
    defs_file.write(" aradius     entry {Aperature radius for photometry}           "   +  aradius + "\n" )
    defs_file.write(" cradius     entry {PSF Clipping radius for optimal photom}    "   +  cradius + "\n" )    
    defs_file.write(" starnumber  entry {Star number to optimize photometry for}    "   +  starnumber + "\n" )  
    defs_file.write(" skyskew     entry {Sky profile skew for optimal photometry}   "   +  skyskew   + "\n" )  
    defs_file.write(" skyfit      entry {Sky fitting switch for optimal photometry} "   +  skyfit    + "\n" )
    defs_file.write(" gain        entry {Instrument Gain e-/ADU}         " +  gain       + "\n" )    
    defs_file.write(" verbose     checkbox {Verbose output switch}       " +  verbose    + "\n" )    
    defs_file.write("endparam\n")

    defs_file.write("\n")
    defs_file.write("param gifparams\n")
    defs_file.write(" wildcard entry {Commong wildcard string for PNG images} " + wildcard + "\n" )
    defs_file.write(" nframes  entry {Number of last N frames to loop} " + nframes + "\n" )
    defs_file.write(" tsleep   entry {Sleep time before checking for new frames [s]} " + tsleep + "\n" )
    defs_file.write("endparam\n")
    defs_file.write("\n")
    defs_file.write("RTPhoS Settings\n")
    defs_file.write("*\n")
    defs_file.write("menu\n")
    defs_file.write("$param(rtphosdef);  $param(gifparams); xterm -hold -sb -sl 2000 -e bash -c \"$rtphos/write_rtphos_defaults.py $rtphos $pathdefs $pathlog $pathans $pathdata $pathbias $pathdark $pathflat $stringbias $stringdark $stringflat $namebias $namedark $nameflat $cprefix $sradius $aradius $cradius $starnumber $skyskew $skyfit $gain $verbose $wildcard $nframes $tsleep \" \n")
    defs_file.write("\n")
    defs_file.write("RUN RTPHOS\n")
    defs_file.write("*\n")
    defs_file.write("menu\n")
    defs_file.write("$param(rtphosdef); xterm -hold -sb -sl 2000 -e \"tcsh -c \\\"python -u $rtphos/run_rtphos.py $xpa_method |& tee $pathlog\\\" \" \n")
    defs_file.write("\n")
    defs_file.write("GIF MONITOR\n")
    defs_file.write("*\n")
    defs_file.write("menu\n")
    defs_file.write("$param(gifparams); $param(rtphosdef); xterm -hold -sb -sl 2000 -e bash -c \" $rtphos/loop_gif.py $filename $wildcard $nframes $tsleep\"\n")

    defs_file.close()

    print "Wrote :"+pathans


if  __name__ == "__main__":

    import sys
    rtphos        = sys.argv[1]
    pathdefs      = sys.argv[2]
    pathlog       = sys.argv[3]
    pathans        = sys.argv[4]
    pathdata       = sys.argv[5]
    pathbias       = sys.argv[6]
    pathdark       = sys.argv[7]
    pathflat       = sys.argv[8]
    stringbias     = sys.argv[9]
    stringdark     = sys.argv[10]
    stringflat     = sys.argv[11]
    namebias       = sys.argv[12]
    namedark       = sys.argv[13]
    nameflat       = sys.argv[14]
    cprefix        = sys.argv[15]
    sradius        = sys.argv[16]
    aradius        = sys.argv[17]
    cradius        = sys.argv[18]
    starnumber     = sys.argv[19]
    skyskew        = sys.argv[20]
    skyfit         = sys.argv[21]
    gain           = sys.argv[22]
    verbose        = sys.argv[23]
    wildcard       = sys.argv[24]
    nframes        = sys.argv[25]
    tsleep         = sys.argv[26]

    write_rtphos_defaults( rtphos, pathdefs, pathlog, pathans, pathdata, pathbias, pathdark, pathflat, stringbias, stringdark, stringflat, namebias, namedark, nameflat, cprefix, sradius, aradius, cradius, starnumber, skyskew, skyfit, gain, verbose, wildcard, nframes, tsleep)
