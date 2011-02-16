#ifndef READ_FITS_H_GUARD
#define READ_FITS_H_GUARD
int nelem(char *filename, long *length);
int read_data(char *filename, char *ch, long firstelem, long nelements, double *data);
#endif
