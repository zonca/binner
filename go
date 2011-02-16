rm *.fits
make
mpirun -n 4 ./binner
