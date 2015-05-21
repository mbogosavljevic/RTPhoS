g95 -c opt_extr.f90 
g95 -o myopphot myopphot.f90 marq.o opt_extr.o subs.o -L/home/zac/Software/cfitsio -lcfitsio

