#include "vectors.h"
#include "fftw3.h"

float *gen_fvect(long length)
{
   float *v;

   v = (float *) fftwf_malloc((size_t) (sizeof(float) * length));
   if (!v) {
      perror("\nError in gen_fvect()");
      printf("\n");
      exit(-1);
   }
   return v;
}


double *gen_dvect(long length)
{
   double *v;

   v = (double *) fftwf_malloc((size_t) (sizeof(double) * length));
   if (!v) {
      perror("\nError in gen_dvect()");
      printf("\n");
      exit(-1);
   }
   return v;
}


short *gen_svect(long length)
{
   short *v;

   v = (short *) fftwf_malloc((size_t) (sizeof(short) * length));
   if (!v) {
      perror("\nError in gen_svect()");
      printf("\n");
      exit(-1);
   }
   return v;
}


int *gen_ivect(long length)
{
   int *v;

   v = (int *) fftwf_malloc((size_t) (sizeof(int) * length));
   if (!v) {
      perror("\nError in gen_ivect()");
      printf("\n");
      exit(-1);
   }
   return v;
}


long *gen_lvect(long length)
{
   long *v;

   v = (long *) fftwf_malloc((size_t) (sizeof(long) * length));
   if (!v) {
      perror("\nError in gen_lvect()");
      printf("\n");
      exit(-1);
   }
   return v;
}


unsigned char *gen_bvect(long length)
{
   unsigned char *v;

   v = (unsigned char *) fftwf_malloc((size_t) (sizeof(unsigned char) * length));
   if (!v) {
      perror("\nError in gen_bvect()");
      printf("\n");
      exit(-1);
   }
   return v;
}
