f95 -c marq.f90
f95 -c opt_extr.f90
f95 -o optimal optimal.f90 marq.o opt_extr.o -L/home/milan/RTPhoS_dependencies/cfitsio/lib -lcfitsio
rm marq.o opt_extr.o

