#include "psrfits.h"
#include "wapp_key.h"

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

long long get_WAPP_info(FILE *files[], int numfiles, int numwapps,
                        struct HEADERP **h, struct wappinfo *w);

void fill_psrfits_struct(int numwapps, struct HEADERP *h, 
                         struct wappinfo *w, struct psrfits *pf);
