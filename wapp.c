#include <stdio.h>
#include <math.h>
#include <string.h>
#include "chkio.h"
#include "vectors.h"
#include "wapp.h"

double slaCldj(int iy, int im, int id, int *j);
static double inv_cerf(double input);
static void vanvleck3lev(float *rho, int npts);
static void vanvleck9lev(float *rho, int npts);

char *get_hdr_string(struct HEADERP *h, char *name, int *slen)
{
    struct HEADERVAL val;
    char *cptr;

    if(find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    cptr = (char *)val.value;
    *slen = strlen(cptr);
    return cptr;
}

double get_hdr_double(struct HEADERP *h, char *name)
{
    struct HEADERVAL val;
    double dval;

    if(find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    dval = *((double *)val.value);
    return dval;
}

int get_hdr_int(struct HEADERP *h, char *name)
{
    struct HEADERVAL val;
    int ival;

    if(find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    ival = *((int *)val.value);
    return ival;
}

long long get_hdr_longlong(struct HEADERP *h, char *name)
{
    struct HEADERVAL val;
    long long llval;

    if(find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    llval = *((long long *)val.value);
    return llval;
}

double *get_hdr_double_arr(struct HEADERP *h, char *name, int *len)
{
    struct HEADERVAL val;
    int ii;
    double *dptr, *darr;

    if(find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    dptr = (double *)val.value;
    *len = val.key->alen;
    /* Note that this needs to be freed! */
    darr = gen_dvect(*len);
    for (ii=0; ii<*len; ii++, dptr++) {
        darr[ii] = *dptr;
    }
    return darr;
}

int *get_hdr_int_arr(struct HEADERP *h, char *name, int *len)
{
    struct HEADERVAL val;
    int ii, *iptr, *iarr;

    if(find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    iptr = (int *)val.value;
    *len = val.key->alen;
    /* Note that this needs to be freed! */
    iarr = gen_ivect(*len);
    for (ii=0; ii<*len; ii++, iptr++) {
        iarr[ii] = *iptr;
    }
    return iarr;
}

static long double UT_strings_to_MJD(char *obs_date, char *start_time)
{
   int year, month, day, hour, min, sec, err;
   long double LMJD;
   
   sscanf(obs_date, "%4d%2d%2d", &year, &month, &day);
   sscanf(start_time, "%2d:%2d:%2d", &hour, &min, &sec);
   LMJD = slaCldj(year, month, day, &err);
   LMJD += (hour + (min + (sec / 60.0)) / 60.0) / 24.0;
   return LMJD;
}

static void set_wappinfo(struct HEADERP *h, struct wappinfo *w)
{
    int ival;
    
    // Header version
    w->header_version = get_hdr_int(h, "header_version");
    
    // Header Length
    w->header_size = h->headlen + h->offset;
    
    // Number of channels (i.e. lags)
    w->numchans = get_hdr_int(h, "num_lags");
    
    // Number of IFs (i.e. polns) in the data
    w->numifs = get_hdr_int(h, "nifs");
    
    // 16- or 32-bit lags
    ival = get_hdr_int(h, "lagformat");
    if (ival == 0)
        w->bits_per_lag = 16;
    else if (ival == 1)
        w->bits_per_lag = 32;
    else
        printf("\nERROR:  Unrecognized number of bits per sample ('lagformat')!\n\n");
    
    // How many bytes are there in a sample in time (i.e. all the lags)
    w->bytes_per_sample = (w->bits_per_lag/8) * w->numifs * w->numchans;

    // 3- or 9-level correlations
    ival = get_hdr_int(h, "level");
    if (ival==1)
        w->corr_level = 3;
    else if (ival==2)
        w->corr_level = 9;
    else
        printf("\nERROR:  Unrecognized 'level' setting!\n\n");
    
    // Is the band inverted?  (i.e. lower sideband)
    w->invertband = 0;
    if (get_hdr_int(h, "freqinversion")) {
        w->invertband = (w->invertband==1 ? 0 : 1);
    }
    if (get_hdr_int(h, "iflo_flip")) {
        w->invertband = (w->invertband==1 ? 0 : 1);
    }
    printf("freqinversion = %d   iflo_flip = %d : using w->invertband = %d\n", 
           get_hdr_int(h, "freqinversion"),
           get_hdr_int(h, "iflo_flip"),
           w->invertband);
    
    // Sampling time in sec
    w->dt = get_hdr_double(h, "wapp_time") / 1.0e6;
    
    // Bandwidth in MHz
    w->BW = get_hdr_double(h, "bandwidth");
    
    // Center freq of the band
    w->fctr = get_hdr_double(h, "cent_freq");

    // Channel BW in MHz
    w->df = fabs(w->BW / w->numchans);

    // Center freq of the lowest freq channel
    w->lofreq = w->fctr - 0.5 * w->BW + 0.5 * w->df;

    // Correlator scaling
    {
        double dtus, corr_rate;
        
        dtus = w->dt * 1.0e6;
        corr_rate = 1.0 / (dtus - WAPP_DEADTIME);
        w->corr_scale = corr_rate / w->BW;
        /* Correction for narrow band use */
        if (w->BW < 50.0)
            w->corr_scale = corr_rate / 50.0;
        if (w->corr_level == 9)      /* 9-level sampling */
            w->corr_scale /= 16.0;
        if (get_hdr_int(h, "sum")) /* summed IFs (search mode) */
            w->corr_scale /= 2.0;
        w->corr_scale *= pow(2.0, (double) get_hdr_int(h, "lagtrunc"));
    }

    // RA (J2000) in decimal degrees
    {
        double hr, m, s, dval;
        dval = get_hdr_double(h, "src_ra");
        hr = (int) floor(dval / 10000.0);
        m = (int) floor((dval - hr * 10000) / 100.0);
        s = dval - hr * 10000 - m * 100;
        w->ra = (hr + (m + s / 60.0) / 60.0) * 15.0;
    }

    // DEC (J2000) in decimal degrees
    {
        double d, m, s, dval;
        dval = get_hdr_double(h, "src_dec");
        d = (int) floor(fabs(dval) / 10000.0);
        m = (int) floor((fabs(dval) - d * 10000) / 100.0);
        s = fabs(dval) - d * 10000 - m * 100;
        w->dec = (d + (m + s / 60.0) / 60.0);
        if (dval < 0.0)
            w->dec = -w->dec;
    }

    // Epoch (MJD) of the first sample
    {
        int len;
        char *cptr1, *cptr2;
        cptr1 = get_hdr_string(h, "obs_date", &len);
        cptr2 = get_hdr_string(h, "start_time", &len);
        w->MJD_epoch = UT_strings_to_MJD(cptr1, cptr2);
    }


#if 1
    printf("header_version = %d\n", w->header_version);
    printf("header_size = %d\n", w->header_size);
    printf("numchans = %d\n", w->numchans);
    printf("numifs = %d\n", w->numifs);
    printf("bits_per_lag = %d\n", w->bits_per_lag);
    printf("corr_level = %d\n", w->corr_level);
    printf("invertband = %d\n", w->invertband);
    printf("dt = %.10f\n", w->dt);
    printf("BW = %f\n", w->BW);
    printf("fctr = %f\n", w->fctr);
    printf("df = %f\n", w->df);
    printf("lofreq = %15.10f\n", w->lofreq);
    printf("corr_scale = %f\n", w->corr_scale);
    printf("ra = %15.10f\n", w->ra);
    printf("dec = %15.10f\n", w->dec);
    printf("MJD_epoch = %Lf\n", w->MJD_epoch);
#endif
}

long long get_WAPP_info(FILE *files[], int numfiles, int numwapps,
                        struct HEADERP *h, struct wappinfo *w)
{
    int ii, jj;
    struct HEADERP *h2;
    struct wappinfo w2;
    long long N=0;
    
    /* Read the header of the first file with the yacc/lex generated tools */
    h = head_parse(files[0]);
    set_wappinfo(h, w);
    close_parse(h);
    
    /* Skip the ASCII and binary headers of the first set of WAPP files */
    for (ii = 0; ii < numwapps; ii++)
        chkfseek(files[ii], w->header_size, SEEK_SET);
    
    // Number of samples in the file
    N = (chkfilelen(files[0], 1) - w->header_size) / w->bytes_per_sample;

    // Now check the other files
    for (ii = 1; ii < numfiles / numwapps; ii++) {

        /* Read the header with the yacc/lex generated tools */
        h2 = head_parse(files[ii*numwapps]);
        set_wappinfo(h, &w2);

        /* Skip the ASCII and binary headers of the other WAPP files */
        for (jj = 0; jj < numwapps; jj++)
            chkfseek(files[ii*numwapps+jj], w->header_size, SEEK_SET);
        close_parse(h2);

        // Do a basic check to see if the files are similar
        if (w2.numchans != w->numchans) {
            printf("Number of channels (file %d) is not the same!\n\n", ii + 1);
        }
        if (w2.numifs != w->numifs) {
            printf("Number of IFs (file %d) is not the same!\n\n", ii + 1);
        }
        if (w2.dt != w->dt) {
            printf("Sample time (file %d) is not the same!\n\n", ii + 1);
        }
        if (w2.df != w->df) {
            printf("Channel width (file %d) is not the same!\n\n", ii + 1);
        }
        
        // Number of samples in the file
        N += (chkfilelen(files[ii*numwapps], 1) - 
              w->header_size) / w->bytes_per_sample;
    }
    return N;
}

