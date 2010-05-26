#include "psrfits.h"
#include "wapp_key.h"
#include "fftw3.h"

/* time between correlator dumps in us */
#define WAPP_DEADTIME 0.34

struct wappinfo {
    int header_version;
    int header_size;
    int numchans;
    int numifs;
    int numsamples;
    int bits_per_lag;
    int bytes_per_sample;
    int corr_level;
    int invertband;
    int beamnum;
    double dt;
    double BW;
    double fctr;
    double df;
    double lofreq;
    double corr_scale;
    double ra;
    double dec;
    long double MJD_epoch;
    char date_obs[80];
};

long long get_WAPP_info(char *filename, FILE *files[], 
                        int numfiles, int numwapps,
                        struct HEADERP **h, struct wappinfo *w);

void fill_psrfits_struct(int numwapps, int numbits, struct HEADERP *h, 
                         struct wappinfo *w, struct psrfits *pf);

int read_WAPP_lags(FILE *infiles[], int numfiles, int numwapps,
                   unsigned char *data, struct wappinfo *w);

void WAPP_lags_to_spectra(int numwapps, struct wappinfo *w, 
                          void *rawdata, float *spectra, float *lags, 
                          fftwf_plan fftplan);

