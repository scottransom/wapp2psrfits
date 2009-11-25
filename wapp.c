#include <stdio.h>
#include <math.h>
#include <string.h>
#include "chkio.h"
#include "vectors.h"
#include "wapp.h"
#include "psrfits.h"
#include "slalib.h"

#ifndef DEGTORAD
#define DEGTORAD 0.017453292519943295769236907684886127134428718885417
#endif
#ifndef RADTODEG
#define RADTODEG 57.29577951308232087679815481410517033240547246656
#endif
#ifndef SOL
#define SOL 299792458.0
#endif

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


static long double UT_strings_to_MJD(char *obs_date, char *start_time, 
                                     char *date_obs)
{
   int year, month, day, hour, min, sec, err;
   double dMJD;
   long double LMJD;
   
   sscanf(obs_date, "%4d%2d%2d", &year, &month, &day);
   sscanf(start_time, "%2d:%2d:%2d", &hour, &min, &sec);
   sprintf(date_obs, "%04d-%02d-%02dT%02d:%02d:%06.3f",
           year, month, day, hour, min, (double) sec);
   slaCaldj(year, month, day, &dMJD, &err);
   LMJD = dMJD + (hour + (min + (sec / 60.0)) / 60.0) / 24.0;
   return LMJD;
}


// Return the beam FWHM in degrees for obs_freq in MHz
// and dish_diam in m
static double beam_FWHM(double obs_freq, double dish_diam)
{
    double lambda = SOL/(obs_freq*1e6);
    return 1.2 * lambda / dish_diam * RADTODEG;
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
        w->MJD_epoch = UT_strings_to_MJD(cptr1, cptr2, w->date_obs);
    }


#if 1
    printf("header_version = %d\n", w->header_version);
    printf("header_size = %d\n", w->header_size);
    printf("numchans = %d\n", w->numchans);
    printf("numifs = %d\n", w->numifs);
    printf("bits_per_lag = %d\n", w->bits_per_lag);
    printf("corr_level = %d\n", w->corr_level);
    printf("invertband = %d\n", w->invertband);
    printf("dt = %.10g\n", w->dt);
    printf("BW = %.10g\n", w->BW);
    printf("fctr = %.10g\n", w->fctr);
    printf("df = %.10g\n", w->df);
    printf("lofreq = %.10g\n", w->lofreq);
    printf("corr_scale = %.10g\n", w->corr_scale);
    printf("ra = %.10g\n", w->ra);
    printf("dec = %.10g\n", w->dec);
    printf("MJD_epoch = %20.15Lf\n", w->MJD_epoch);
    printf("date_obs = '%s'\n", w->date_obs);
#endif
}


static int compare_wapp_files_basic(int filenum, struct wappinfo *w1, 
                                    struct wappinfo *w2)
{
    int good=1;
    
    if (w2->numchans != w1->numchans) {
        printf("Error:  Number of channels (%d vs %d) differs between files %d and 1!\n\n", 
               w2->numchans, w1->numchans, filenum+1);
        good = 0;
    }
    if (w2->numifs != w1->numifs) {
        printf("Error:  Number of IFs (%d vs %d) differs between files %d and 1!\n\n", 
               w2->numifs, w1->numifs, filenum+1);
        good = 0;
    }
    if (w2->dt != w1->dt) {
        printf("Error:  Sample time (%.4f vs %.4f) differs between files %d and 1!\n\n", 
               w2->dt, w1->dt, filenum+1);
        good = 0;
    }
    if (w2->df != w1->df) {
        printf("Error:  Channel width (%.4f vs %.4f) differs between files %d and 1!\n\n", 
               w2->df, w1->df, filenum+1);
        good = 0;
    }
    return good;
}

static int compare_samewapp_files(int filenum, struct wappinfo *w1, 
                                  struct wappinfo *w2)
{
    int good = compare_wapp_files_basic(filenum, w1, w2);
    
    if (w2->fctr != w1->fctr) {
        printf("Error:  Center freq (%.4f vs %.4f) differs between files %d and 1!\n\n", 
               w2->fctr, w1->fctr, filenum+1);
        good = 0;
    }
    return good;
}
        

static int compare_diffwapp_files(int filenum, int numwapps, int wappindex, 
                                  struct wappinfo *w1, struct wappinfo *w2)
{
    int good = compare_wapp_files_basic(filenum, w1, w2);
    double diff_BWs;
    
    if (filenum < numwapps) {
        // Compare the starting epochs
        if (fabs(w2->MJD_epoch - w1->MJD_epoch) * 86400.0 > 1.e-5) {
            printf("Error:  Epoch (%20.15Lf vs %20.15Lf) differs between files %d and 1!\n\n", 
                   w2->MJD_epoch, w1->MJD_epoch, filenum+1);
            good = 0;
        }
    }
    // Compare the center frequencies of the files
    diff_BWs = fabs(w1->fctr - w2->fctr) / w1->BW;
    if (fabs(diff_BWs - wappindex) > 1e-3) {
        printf("Error:  Center freqs (%.2f MHz vs %.2f MHz) are not separated\n"
               "        by the BW (%.2f MHz), in files %d and 1!\n\n", 
               w2->fctr, w1->fctr, w1->BW, filenum+1);
        good = 0;
    }
    return good;
}


static void dec2hms(char *out, double in, int sflag) {
    int sign = 1;
    char *ptr = out;
    int h, m;
    double s;
    if (in<0.0) { sign = -1; in = fabs(in); }
    h = (int)in; in -= (double)h; in *= 60.0;
    m = (int)in; in -= (double)m; in *= 60.0;
    s = in;
    if (sign==1 && sflag) { *ptr='+'; ptr++; }
    else if (sign==-1) { *ptr='-'; ptr++; }
    sprintf(ptr, "%2.2d:%2.2d:%07.4f", h, m, s);
}


void fill_psrfits_struct(int numwapps, struct HEADERP *h, 
                         struct wappinfo *w, struct psrfits *pf)
{
    int slen, ii;
    char *cptr;

    pf->filenum = 0;  // Crucial for initialization
    pf->hdr.nsblk = (int) (1.0/w->dt); // _might_ be a problem...

    // Now set values for our hdrinfo structure
    strcpy(pf->hdr.telescope, "Arecibo");
    cptr = get_hdr_string(h, "obs_type", &slen);
    if (strncmp("PULSAR_SEARCH", cptr, slen)==0) {
        strcpy(pf->hdr.obs_mode, "SEARCH");
    } else {
        printf("Error:  Wapp data is not in search format!\n\n");
        exit(1);
    }
    strcpy(pf->hdr.backend, "WAPP");
    cptr = get_hdr_string(h, "frontend", &slen);
    strncpy(pf->hdr.frontend, cptr, slen);
    cptr = get_hdr_string(h, "observers", &slen);
    strncpy(pf->hdr.observer, cptr, slen);
    cptr = get_hdr_string(h, "project_id", &slen);
    strncpy(pf->hdr.project_id, cptr, slen);
    cptr = get_hdr_string(h, "src_name", &slen);
    strncpy(pf->hdr.source, cptr, slen);
    strcpy(pf->hdr.date_obs, w->date_obs);
    pf->hdr.scanlen = get_hdr_double(h, "obs_time");

    strcpy(pf->hdr.poln_type, "LIN"); // set based on known receivers
    if (get_hdr_int(h, "sum")) {
        strcpy(pf->hdr.poln_order, "AA+BB");
        pf->hdr.summed_polns = 1;
    } else if (w->numifs==1) {
        strcpy(pf->hdr.poln_order, "AA");
        pf->hdr.summed_polns = 0;
    }
    strcpy(pf->hdr.track_mode, "TRACK");  // Potentially not-true?
    strcpy(pf->hdr.cal_mode, "OFF");      // Potentially not-true?
    strcpy(pf->hdr.feed_mode, "FA"); // check this...

    if (get_hdr_int(h, "isalfa")) {
        // TODO:
        //   Need to set the beam number here...
        //   Also should correct positions and paralactic angles etc...
    }
 
    pf->hdr.dt = w->dt;
    pf->hdr.fctr = w->fctr;
    pf->hdr.BW = w->BW;
    pf->hdr.beam_FWHM = beam_FWHM(pf->hdr.fctr, 300.0);
    pf->hdr.nchan = w->numchans;
    pf->hdr.orig_nchan = w->numchans;
    pf->hdr.orig_df = pf->hdr.df = pf->hdr.BW / pf->hdr.nchan;
    pf->hdr.nbits = 8;  // This needs to be re-set by the command line option
    pf->hdr.npol = w->numifs;
    pf->hdr.MJD_epoch = w->MJD_epoch;
    pf->hdr.start_day = (int) (w->MJD_epoch);
    pf->hdr.start_sec = (w->MJD_epoch - pf->hdr.start_day) * 86400.0;
    pf->hdr.scan_number = get_hdr_int(h, "scan_number");
    pf->hdr.ra2000 = w->ra;
    dec2hms(pf->hdr.ra_str, pf->hdr.ra2000/15.0, 0);
    pf->hdr.dec2000 = w->dec;
    dec2hms(pf->hdr.dec_str, pf->hdr.dec2000, 1);
    pf->hdr.azimuth = get_hdr_double(h, "start_az");
    pf->hdr.zenith_ang = get_hdr_double(h, "start_za");
    pf->hdr.rcvr_polns = 2;
    pf->hdr.offset_subint = 0;
    pf->hdr.onlyI = 0;
    pf->hdr.ds_time_fact = 1;
    pf->hdr.ds_freq_fact = 1;
    pf->hdr.chan_dm = 0.0;
    pf->hdr.fd_hand = pf->hdr.be_phase = 0;  // This is almost certainly not correct
    pf->hdr.fd_sang = pf->hdr.fd_xyph = 0.0; // This is almost certainly not correct
    pf->hdr.feed_angle = 0.0;            // This is almost certainly not correct
    pf->hdr.cal_freq = pf->hdr.cal_dcyc = pf->hdr.cal_phs = 0.0;  //  ditto

    // Now set values for our subint structure
    pf->sub.tel_az = get_hdr_double(h, "start_az");
    pf->sub.tel_zen = get_hdr_double(h, "start_za");
    pf->sub.lst = get_hdr_double(h, "start_lst");
    pf->sub.tsubint = pf->hdr.nsblk * pf->hdr.dt;
    pf->sub.ra = pf->hdr.ra2000;
    pf->sub.dec = pf->hdr.dec2000;
    pf->sub.offs = 0.5 * pf->sub.tsubint;
    slaEqgal(pf->hdr.ra2000*DEGTORAD, pf->hdr.dec2000*DEGTORAD, 
             &pf->sub.glon, &pf->sub.glat);
    pf->sub.glon *= RADTODEG;
    pf->sub.glat *= RADTODEG;
    pf->sub.feed_ang = 0.0;  // These are wrong too...
    pf->sub.pos_ang = 0.0;   // These are wrong too...
    pf->sub.par_ang = 0.0;   // These are wrong too...
    pf->sub.bytes_per_subint = (pf->hdr.nbits * pf->hdr.nchan * 
                               pf->hdr.npol * pf->hdr.nsblk) / 8;
    pf->sub.FITS_typecode = TBYTE;  // 11 = byte

    // Create and initialize the subint arrays
    pf->sub.dat_freqs = gen_fvect(pf->hdr.nchan * numwapps);
    pf->sub.dat_weights = gen_fvect(pf->hdr.nchan * numwapps);
    for (ii = 0 ; ii < pf->hdr.nchan * numwapps; ii++) {
        pf->sub.dat_freqs[ii] = w->lofreq + ii * pf->hdr.df;
        pf->sub.dat_weights[ii] = 1.0;
    }

    // The following need to be adjusted for 4-bit data at least...
    pf->sub.dat_offsets = gen_fvect(pf->hdr.nchan * pf->hdr.npol * numwapps);
    pf->sub.dat_scales = gen_fvect(pf->hdr.nchan * pf->hdr.npol * numwapps);
    for (ii = 0 ; ii < pf->hdr.nchan * pf->hdr.npol * numwapps; ii++) {
        pf->sub.dat_offsets[ii] = 0.0;
        pf->sub.dat_scales[ii] = 1.0;
    }
 
    // This is the raw data block that will be updated 
    // for each row of the PSRFITS file
    pf->sub.data = gen_bvect(pf->sub.bytes_per_subint);
}


long long get_WAPP_info(FILE *files[], int numfiles, int numwapps,
                        struct HEADERP **h, struct wappinfo *w)
{
    int ii, wappindex;
    struct HEADERP *h2;
    struct wappinfo w2;
    long long N=0;
    
    // Read the header of the first file with the yacc/lex generated tools
    // This sets the basic parameters of the conversion
    *h = head_parse(files[0]);
    set_wappinfo(*h, w);
    // Number of samples in the file
    w->numsamples = (chkfilelen(files[0], 1) - w->header_size) / 
        w->bytes_per_sample;
    // This will be the total number of samples to convert
    N = w->numsamples;

    // Skip the ASCII and binary header
    chkfseek(files[0], w->header_size, SEEK_SET);
    
    // loop through all the other files and check/prep them
    for (ii = 1; ii < numfiles; ii++) {
       
        // Read the header with the yacc/lex generated tools
        h2 = head_parse(files[ii]);
        set_wappinfo(h2, &w2);
        close_parse(h2);
        // Number of samples in the file
        w2.numsamples = (chkfilelen(files[ii], 1) - w2.header_size) / 
            w2.bytes_per_sample;
        
        // Skip the ASCII and binary header
        chkfseek(files[ii], w2.header_size, SEEK_SET);
        
        // Do a basic check to see if the files are similar
        wappindex = ii % numwapps;
        if (wappindex == 0) { // Same WAPP
            if (!compare_samewapp_files(ii, w, &w2)) 
                exit(1);
            N += w2.numsamples;
        } else { // Different WAPPs
            if (!compare_diffwapp_files(ii, numwapps, wappindex, w, &w2)) 
                exit(1);
        }
    }
    return N;
}

