#ifndef READ_FITS_LIB_H_GUARD
#define READ_FITS_LIB_H_GUARD
int latest_exchange(int od, int frequency, char type, const char* path, char* filename);
int nelem(const char* filename, long *length);
int read_data(const char *filename, const char *extname, int colnum, long firstelem, long nelements, double *data);
int write_map(const char *filename, double *intensity, double *q, double *u, long nside);
#endif
