f95 -fcheck=all -c opt_extr.f90 marq.f90
f95 -fcheck=all -o optimal optimal.f90 marq.o opt_extr.o -L/home/zac/Software/cfitsio -lcfitsio
rm marq.o opt_extr.o

