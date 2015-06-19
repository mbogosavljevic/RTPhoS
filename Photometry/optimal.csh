g95 -c opt_extr.f90 marq.f90
g95 -o optimal optimal.f90 marq.o opt_extr.o -L/home/zac/Software/cfitsio -lcfitsio
rm marq.o opt_extr.o subs.o

