#include <string.h>
#include <math.h>
#include <stdio.h>
#include "fitsio.h"
#include <dirent.h>
#include <stdbool.h>


const size_t NUM_P = 10;

void swap( void **p1,  void **p2)
{
  void *pt = *p1;
  *p1 = *p2;
  *p2 = pt;
}

void str_sort(const char *p[], int n)
{
  char *pTemp = NULL;                 /* Temporary pointer               */
  bool sorted = false;                /* Strings sorted indicator        */
  while(!sorted)                      /* Loop until there are no swaps   */
  {
    sorted = true;                    /* Initialize to indicate no swaps */
    int i;
    for(i = 0 ; i<n-1 ; i++ )
      if(strcmp(p[i], p[i + 1]) > 0)
      {
        sorted = false;               /* indicate we are out of order    */
        swap(&p[i], &p[i+1]);         /* Swap the pointers               */
      }
  }
}

int latest_exchange(int od, int frequency, char type, const char* path, char* filename)
{
    DIR *dp;
    struct dirent *ep;

    sprintf(filename, "%s/%04d/",path,od);

    char namestart[20];
    sprintf(namestart, "L%03d-%04d-%c-",frequency,od,type);

    char *pS[NUM_P];                     /* Array of string pointers         */
    
    dp = opendir (filename);
    int i=0;
    if (dp != NULL)
     {
       while (ep = readdir (dp))
           if (strncmp(ep->d_name, namestart,12) == 0) {
         puts (ep->d_name);
            pS[i] = (char*)malloc(strlen( ep->d_name) + 1);
            strcpy(pS[i], ep->d_name); 
            /* DEBUG
            printf("%d\n",i);
            */
            i++;
           }
       (void) closedir (dp);
     }
    else
     perror ("Couldn't open the directory");
    str_sort( pS, i);                          /* Sort strings           */
    /* DEBUG
    int j=0;
    for (j=0 ; j<i;j++)
    {
        printf("%d\n",j);
        puts(pS[j]);
    }
    */

    strcat(filename, pS[i-1]);
    return 0;
}
    

int nelem(const char* filename, long *length)
{
    int status = 0;
    fitsfile *fptr;         
    fits_open_file(&fptr, filename, READONLY, &status);
    fprintf(stderr, "%s :reading.\n", filename);
    fits_movnam_hdu(fptr, BINARY_TBL, "DATA", 0, &status);
    fits_read_key(fptr, TLONG, "NAXIS2", length , NULL, &status);
    fits_close_file(fptr, &status);
    if (status)          /* print any error messages */
        fits_report_error(stderr, status);
    return(status);
}

int read_data(const char *filename,const char *extname, int colnum, long firstelem, long nelements, double *data)
{
    firstelem++; //from 0 to 1 based
    int status = 0, anynul;
    double doublenull = 0.;
    fitsfile *fptr;         
    fits_open_file(&fptr, filename, READONLY, &status);
    fits_movnam_hdu(fptr, BINARY_TBL, extname, 0, &status);
    fits_read_col(fptr, TDOUBLE, colnum, firstelem, 1, nelements, &doublenull, data, &anynul, &status);
    fits_close_file(fptr, &status);
    if (status)          /* print any error messages */
        fits_report_error(stderr, status);
    return(status);
}

int write_map(const char *filename, double *q, double *u, long nrows)
{
    
    fitsfile *fptr;
    int status = 0;
    int hdutype;
    int tfields=2;
    char extname[] = "BINTABLE";   /* extension name */
    char *ttype[] = { "Q","U" };
    char *tform[] = { "1E","1E" };
    char *tunit[] = { " ", " " };
    if (fits_create_file(&fptr, filename, &status)) 
    fprintf(stderr, "%s (%d): Could not create new fits file.\n", 
        __FILE__, __LINE__);

    //double first[maplength];
    //double second[maplength];

    //printf("%d\n",maplength);
    //int i = 0;
    //for (i = 0; i < nrows;i=i+2){
    //    first[i/2] = signal[i];
    //    second[i/2 + 1] = signal[i +1];
    //}
    /* move to 1st HDU  */
    fits_movabs_hdu(fptr, 1, &hdutype, &status);
    //long nrows = 12L*nside*nside*2;
    long nside = sqrtl(nrows/12L);

    /* append a new empty binary table onto the FITS file */
    fits_create_tbl( fptr, BINARY_TBL, nrows, tfields, ttype, tform,
            tunit, extname, &status);

    fits_write_key(fptr, TSTRING, "ORDERING", "RING", 
             "Pixel ordering scheme, either RING or NESTED", &status);

    fits_write_key(fptr, TLONG, "NSIDE", &nside, "Resolution parameter for HEALPIX", &status);
    long firstrow  = 1;  /* first row in table to write   */
    long firstelem = 1;  /* first element in row  (ignored in ASCII tables)  */

    if (fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, nrows, q,
             &status))
    fprintf(stderr, "%s (%d): Could not write signal.\n", __FILE__, __LINE__);
    if (fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, nrows, u,
             &status))
    fprintf(stderr, "%s (%d): Could not write signal.\n", __FILE__, __LINE__);

    fits_close_file(fptr, &status);

    return true;
}
