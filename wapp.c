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
#ifndef SWAP
/* Swaps two variables of undetermined type */
#define SWAP(a,b) tempzz=(a);(a)=(b);(b)=tempzz;
#endif

static double inv_cerf(double input);
static void vanvleck3lev(float *rho, int npts);
static void vanvleck9lev(float *rho, int npts);

int wapp_beamnum(char *str)
{
    int ii = -1;
    // Check if filename ends in "_N.wapp"
    // Return -1 if the file is not a WAPP file or wrong beam num
    if ((0 == strncmp(".wapp", str + strlen(str) - 5, 5)) &&
        (str[strlen(str) - 7] == '_')) {
        ii = atoi(str + strlen(str) - 6);
        if (ii > 7)
            ii = -1;            // Only 8 WAPPs...
    }
    return ii;
}

char *get_hdr_string(struct HEADERP *h, char *name, int *slen)
{
    struct HEADERVAL val;
    char *cptr;

    if (find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    cptr = (char *) val.value;
    // Is the following + 1 really correct?  I think it causes it
    // to always include the null terminator... SMR
    *slen = strlen(cptr) + 1;
    // printf("'%s' = '%s' (%d chars)\n", name, cptr, *slen);
    return cptr;
}

double get_hdr_double(struct HEADERP *h, char *name)
{
    struct HEADERVAL val;
    double dval;

    if (find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    dval = *((double *) val.value);
    return dval;
}

int get_hdr_int(struct HEADERP *h, char *name)
{
    struct HEADERVAL val;
    int ival;

    if (find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    ival = *((int *) val.value);
    return ival;
}

long long get_hdr_longlong(struct HEADERP *h, char *name)
{
    struct HEADERVAL val;
    long long llval;

    if (find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    llval = *((long long *) val.value);
    return llval;
}

double *get_hdr_double_arr(struct HEADERP *h, char *name, int *len)
{
    struct HEADERVAL val;
    int ii;
    double *dptr, *darr;

    if (find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    dptr = (double *) val.value;
    *len = val.key->alen;
    /* Note that this needs to be freed! */
    darr = gen_dvect(*len);
    for (ii = 0; ii < *len; ii++, dptr++) {
        darr[ii] = *dptr;
    }
    return darr;
}

int *get_hdr_int_arr(struct HEADERP *h, char *name, int *len)
{
    struct HEADERVAL val;
    int ii, *iptr, *iarr;

    if (find_hdrval(h, name, &val)) {
        printf("ERROR:  Can't find '%s' in the WAPP header!\n", name);
        exit(0);
    }
    iptr = (int *) val.value;
    *len = val.key->alen;
    /* Note that this needs to be freed! */
    iarr = gen_ivect(*len);
    for (ii = 0; ii < *len; ii++, iptr++) {
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
    double lambda = SOL / (obs_freq * 1e6);
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
    if (w->numifs > 1) {
        printf("\nERROR:  wapp2psrfits cannot (yet?) handle more than 1 IF!\n\n");
        exit(1);
    }
    // 16- or 32-bit lags
    ival = get_hdr_int(h, "lagformat");
    if (ival == 0)
        w->bits_per_lag = 16;
    else if (ival == 1)
        w->bits_per_lag = 32;
    else {
        printf("\nERROR:  Unrecognized bits per sample (%d, 'lagformat')!\n\n",
               w->bits_per_lag);
        exit(1);
    }

    // How many bytes are there in a sample in time (i.e. all the lags)
    w->bytes_per_sample = (w->bits_per_lag / 8) * w->numifs * w->numchans;

    // 3- or 9-level correlations
    ival = get_hdr_int(h, "level");
    if (ival == 1)
        w->corr_level = 3;
    else if (ival == 2)
        w->corr_level = 9;
    else {
        printf("\nERROR:  Unrecognized 'level' setting! (%d)\n\n", w->corr_level);
        exit(1);
    }

    // Is the band inverted?  (i.e. lower sideband)
    w->invertband = 0;
    if (get_hdr_int(h, "freqinversion")) {
        w->invertband = (w->invertband == 1 ? 0 : 1);
    }
    // These two parameters are cumulative...
    if (w->header_version >= 7) {
        if (get_hdr_int(h, "iflo_flip")) {
            w->invertband = (w->invertband == 1 ? 0 : 1);
        }
    }

    // Sampling time in sec
    w->dt = get_hdr_double(h, "wapp_time") / 1.0e6;

    // Bandwidth in MHz
    w->BW = get_hdr_double(h, "bandwidth");

    // Center freq of the band
    w->fctr = get_hdr_double(h, "cent_freq");

    // Channel BW in MHz
    w->df = fabs(w->BW / w->numchans);

    // Center freq of the lowest freq channel
    // See:  http://www.cv.nrao.edu/~pdemores/wapp/
    // Note:  If band is inverted, since we reorder the channels
    //        explicitly, we don't need to use the negative "B"
    //        factor from Paul's webpage above.  And the following
    //        is correct for USB or LSB data.
    w->lofreq = w->fctr - 0.5 * w->BW + 0.5 * w->df;

    // Correlator scaling
    {
        double dtus, corr_rate;

        dtus = w->dt * 1.0e6;
        corr_rate = 1.0 / (dtus - WAPP_DEADTIME);
        w->corr_scale = corr_rate / w->BW;
        // Correction for narrow band use
        if (w->BW < 50.0)
            w->corr_scale = corr_rate / 50.0;
        if (w->corr_level == 9)     // 9-level sampling
            w->corr_scale /= 16.0;
        if (get_hdr_int(h, "sum"))  // summed IFs (search mode)
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
}


static int compare_wapp_files_basic(int filenum, struct wappinfo *w1,
                                    struct wappinfo *w2)
{
    int good = 1;

    if (w2->numchans != w1->numchans) {
        printf
            ("Error:  Number of channels (%d vs %d) differs between files %d and 1!\n\n",
             w2->numchans, w1->numchans, filenum + 1);
        good = 0;
    }
    if (w2->numifs != w1->numifs) {
        printf
            ("Error:  Number of IFs (%d vs %d) differs between files %d and 1!\n\n",
             w2->numifs, w1->numifs, filenum + 1);
        good = 0;
    }
    if (w2->dt != w1->dt) {
        printf
            ("Error:  Sample time (%.4f vs %.4f) differs between files %d and 1!\n\n",
             w2->dt, w1->dt, filenum + 1);
        good = 0;
    }
    if (w2->df != w1->df) {
        printf
            ("Error:  Channel width (%.4f vs %.4f) differs between files %d and 1!\n\n",
             w2->df, w1->df, filenum + 1);
        good = 0;
    }
    return good;
}

static int compare_samewapp_files(int filenum, struct wappinfo *w1,
                                  struct wappinfo *w2)
{
    int good = compare_wapp_files_basic(filenum, w1, w2);

    if (w2->fctr != w1->fctr) {
        printf
            ("Error:  Center freq (%.4f vs %.4f) differs between files %d and 1!\n\n",
             w2->fctr, w1->fctr, filenum + 1);
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
            printf
                ("Error:  Epoch (%20.15Lf vs %20.15Lf) differs between files %d and 1!\n\n",
                 w2->MJD_epoch, w1->MJD_epoch, filenum + 1);
            good = 0;
        }
    }
    // Compare the center frequencies of the files
    diff_BWs = fabs(w1->fctr - w2->fctr) / w1->BW;
    if (fabs(diff_BWs - wappindex) > 1e-3) {
        printf("Error:  Center freqs (%.2f MHz vs %.2f MHz) are not separated\n"
               "        by the BW (%.2f MHz), in files %d and 1!\n\n",
               w2->fctr, w1->fctr, w1->BW, filenum + 1);
        good = 0;
    }
    return good;
}


static void dec2hms(char *out, double in, int sflag)
{
    int sign = 1;
    char *ptr = out;
    int h, m;
    double s;
    if (in < 0.0) {
        sign = -1;
        in = fabs(in);
    }
    h = (int) in;
    in -= (double) h;
    in *= 60.0;
    m = (int) in;
    in -= (double) m;
    in *= 60.0;
    s = in;
    if (sign == 1 && sflag) {
        *ptr = '+';
        ptr++;
    } else if (sign == -1) {
        *ptr = '-';
        ptr++;
    }
    sprintf(ptr, "%2.2d:%2.2d:%07.4f", h, m, s);
}


void fill_psrfits_struct(int numwapps, int numbits, struct HEADERP *h,
                         struct wappinfo *w, struct psrfits *pf)
{
    int slen, ii;
    char *cptr;

    pf->filenum = 0;            // Crucial for initialization
    pf->hdr.nsblk = (int) (1.0 / w->dt);        // _might_ be a problem...

    // Now set values for our hdrinfo structure
    strcpy(pf->hdr.telescope, "Arecibo");
    cptr = get_hdr_string(h, "obs_type", &slen);
    if (strncmp("PULSAR_SEARCH", cptr, slen) == 0) {
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

    strcpy(pf->hdr.poln_type, "LIN");   // set based on known receivers
    if (get_hdr_int(h, "sum")) {
        strcpy(pf->hdr.poln_order, "AA+BB");
        pf->hdr.summed_polns = 1;
    } else if (w->numifs == 1) {
        strcpy(pf->hdr.poln_order, "AA");
        pf->hdr.summed_polns = 0;
    }
    strcpy(pf->hdr.track_mode, "TRACK");  // Potentially not-true?
    strcpy(pf->hdr.cal_mode, "OFF");      // Potentially not-true?
    strcpy(pf->hdr.feed_mode, "FA");      // check this...

    pf->hdr.beamnum = 0;
    if (get_hdr_int(h, "isalfa"))
        pf->hdr.beamnum = w->beamnum;

    pf->hdr.dt = w->dt;
    pf->hdr.fctr = w->fctr + 0.5 * (numwapps - 1.0) * w->BW;
    pf->hdr.BW = w->BW * numwapps;
    pf->hdr.beam_FWHM = beam_FWHM(pf->hdr.fctr, 300.0);
    pf->hdr.nchan = w->numchans * numwapps;
    pf->hdr.orig_nchan = w->numchans * numwapps;
    pf->hdr.orig_df = pf->hdr.df = pf->hdr.BW / pf->hdr.nchan;
    pf->hdr.nbits = numbits;
    pf->hdr.npol = w->numifs;
    pf->hdr.MJD_epoch = w->MJD_epoch;
    pf->hdr.start_day = (int) (w->MJD_epoch);
    pf->hdr.start_sec = (w->MJD_epoch - pf->hdr.start_day) * 86400.0;
    pf->hdr.scan_number = get_hdr_int(h, "scan_number");
    pf->hdr.ra2000 = w->ra;
    dec2hms(pf->hdr.ra_str, pf->hdr.ra2000 / 15.0, 0);
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
    pf->hdr.fd_hand = pf->hdr.be_phase = 0;     // This is almost certainly not correct
    pf->hdr.fd_sang = pf->hdr.fd_xyph = 0.0;    // This is almost certainly not correct
    pf->hdr.feed_angle = 0.0;   // This is almost certainly not correct
    pf->hdr.cal_freq = pf->hdr.cal_dcyc = pf->hdr.cal_phs = 0.0;        //  ditto

    // Now set values for our subint structure
    pf->sub.tel_az = get_hdr_double(h, "start_az");
    pf->sub.tel_zen = get_hdr_double(h, "start_za");
    pf->sub.lst = get_hdr_double(h, "start_lst");
    pf->sub.tsubint = pf->hdr.nsblk * pf->hdr.dt;
    pf->sub.ra = pf->hdr.ra2000;
    pf->sub.dec = pf->hdr.dec2000;
    pf->sub.offs = 0.5 * pf->sub.tsubint;
    slaEqgal(pf->hdr.ra2000 * DEGTORAD, pf->hdr.dec2000 * DEGTORAD,
             &pf->sub.glon, &pf->sub.glat);
    pf->sub.glon *= RADTODEG;
    pf->sub.glat *= RADTODEG;
    // The following three are unknown or hard to get, I think (SMR)
    pf->sub.feed_ang = 0.0;
    pf->sub.pos_ang = 0.0;
    pf->sub.par_ang = 0.0;
    pf->sub.bytes_per_subint = (pf->hdr.nbits * pf->hdr.nchan *
                                pf->hdr.npol * pf->hdr.nsblk) / 8;
    pf->sub.FITS_typecode = TBYTE;      // 11 = byte

    // Create and initialize the subint arrays
    pf->sub.dat_freqs = gen_fvect(pf->hdr.nchan);
    pf->sub.dat_weights = gen_fvect(pf->hdr.nchan);
    for (ii = 0; ii < pf->hdr.nchan; ii++) {
        pf->sub.dat_freqs[ii] = w->lofreq + ii * pf->hdr.df;
        pf->sub.dat_weights[ii] = 1.0;
    }

    // The following are re-set to try to preserve the band shape later
    pf->sub.dat_offsets = gen_fvect(pf->hdr.nchan * pf->hdr.npol);
    pf->sub.dat_scales = gen_fvect(pf->hdr.nchan * pf->hdr.npol);
    for (ii = 0; ii < pf->hdr.nchan * pf->hdr.npol; ii++) {
        pf->sub.dat_offsets[ii] = 0.0;
        pf->sub.dat_scales[ii] = 1.0;
    }

    // This is the raw data block that will be updated 
    // for each row of the PSRFITS file
    pf->sub.data = gen_bvect(pf->sub.bytes_per_subint);
}


long long get_WAPP_info(char *filename, FILE * files[], int numfiles, int numwapps,
                        struct HEADERP **h, struct wappinfo *w)
{
    int ii, wappindex;
    struct HEADERP *h2;
    struct wappinfo w2;
    long long N = 0;

    // Read the header of the first file with the yacc/lex generated tools
    // This sets the basic parameters of the conversion
    *h = head_parse(files[0]);
    set_wappinfo(*h, w);
    if (get_hdr_int(*h, "isalfa")) {
        w->beamnum = wapp_beamnum(filename);
        if (w->beamnum == -1) {
            printf("Warning!  isalfa = 1, but beamnum is not set!\n");
            exit(1);
        }
    }
    // Number of samples in the file
    w->numsamples = (chkfilelen(files[0], 1) - w->header_size) / w->bytes_per_sample;
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
        if (wappindex == 0) {   // Same WAPP
            if (!compare_samewapp_files(ii, w, &w2))
                exit(1);
            N += w2.numsamples;
        } else {                // Different WAPPs
            if (!compare_diffwapp_files(ii, numwapps, wappindex, w, &w2))
                exit(1);
        }
    }
    return N;
}


int read_WAPP_lags(FILE * infiles[], int numfiles, int numwapps,
                   unsigned char *data, struct wappinfo *w)
// This routine reads a set of lags from the input files *infiles
// which contain 16 or 32 bit lag data from the WAPP correlator at
// Arecibo.  A set of WAPP lags is bytes_per_sample * numchans * #bits
// * numwapps long.  *data must be at least that long
{
    int ii;
    static int currentfile = 0, numread = 0;

    // Make sure our current file number is valid
    if (currentfile >= numfiles / numwapps)
        return 0;

    // First, attempt to read data from the current file
    if (chkfread(data, w->bytes_per_sample, 1, infiles[currentfile * numwapps])) {      // Got data
        // Get data from other WAPPs
        for (ii = 1; ii < numwapps; ii++)
            chkfread(data + ii * w->bytes_per_sample, w->bytes_per_sample, 1,
                     infiles[currentfile * numwapps + ii]);
        numread++;
        return 1;
    } else {                    // Didn't get data
        if (feof(infiles[currentfile * numwapps])) {    // End of file?
            currentfile++;
            return read_WAPP_lags(infiles, numfiles, numwapps, data, w);
        } else {
            printf("\nProblem reading record from WAPP data file:\n");
            printf("   currentfile = %d.  Exiting.\n", currentfile);
            exit(1);
        }
    }
}

void print_WAPP_hdr(struct HEADERP *hdr)
/* Output a WAPP header in human readable form */
{
    int len, vers, inverted = 0;
    char datetime[100];

    vers = get_hdr_int(hdr, "header_version");
    printf("\n             Header version = %d\n", vers);
    printf("        Header size (bytes) = %d\n", hdr->headlen + hdr->offset);
    printf("                Source Name = %s\n", get_hdr_string(hdr, "src_name", &len));
    printf("           Observation Type = %s\n", get_hdr_string(hdr, "obs_type", &len));
    printf(" Observation Date (YYYMMDD) = %s\n", get_hdr_string(hdr, "obs_date", &len));
    printf("    Obs Start UT (HH:MM:SS) = %s\n", get_hdr_string(hdr, "start_time", &len));
    printf("             MJD start time = %.12Lf\n",
           UT_strings_to_MJD(get_hdr_string(hdr, "obs_date", &len),
                             get_hdr_string(hdr, "start_time", &len),
                             datetime));
    printf("                 Project ID = %s\n", get_hdr_string(hdr, "project_id", &len));
    printf("                  Observers = %s\n", get_hdr_string(hdr, "observers", &len));
    printf("                Scan Number = %d\n", get_hdr_int(hdr, "scan_number"));
    printf("    RA (J2000, HHMMSS.SSSS) = %.4f\n", get_hdr_double(hdr, "src_ra"));
    printf("   DEC (J2000, DDMMSS.SSSS) = %.4f\n", get_hdr_double(hdr, "src_dec"));
    printf("        Start Azimuth (deg) = %-17.15g\n", get_hdr_double(hdr, "start_az"));
    printf("     Start Zenith Ang (deg) = %-17.15g\n", get_hdr_double(hdr, "start_za"));
    printf("            Start AST (sec) = %-17.15g\n", get_hdr_double(hdr, "start_ast"));
    printf("            Start LST (sec) = %-17.15g\n", get_hdr_double(hdr, "start_lst"));
    printf("           Obs Length (sec) = %-17.15g\n", get_hdr_double(hdr, "obs_time"));
    printf("      Requested T_samp (us) = %-17.15g\n", get_hdr_double(hdr, "samp_time"));
    printf("         Actual T_samp (us) = %-17.15g\n", get_hdr_double(hdr, "wapp_time"));
    printf("         Central freq (MHz) = %-17.15g\n", get_hdr_double(hdr, "cent_freq"));
    printf("      Total Bandwidth (MHz) = %-17.15g\n", get_hdr_double(hdr, "bandwidth"));
    printf("             Number of lags = %d\n", get_hdr_int(hdr, "num_lags"));
    printf("              Number of IFs = %d\n", get_hdr_int(hdr, "nifs"));
    printf("    Samples since obs start = %lld\n", get_hdr_longlong(hdr, "timeoff"));
    printf("   Other information:\n");
    if (get_hdr_int(hdr, "sum") == 1) printf("      IFs are summed.\n");
    if (get_hdr_int(hdr, "freqinversion")) {
        inverted = (inverted == 1 ? 0 : 1);
    }
    if (vers >= 7) {
        if (get_hdr_int(hdr, "iflo_flip")) {
            inverted = (inverted == 1 ? 0 : 1);
        }
    }
    if (inverted) printf("      Frequency band is inverted.\n");
    if (get_hdr_int(hdr, "lagformat") == 0)
        printf("      Lags are 16 bit integers.\n");
    else
        printf("      Lags are 32 bit integers.\n");
}


void WAPP_lags_to_spectra(int numwapps, struct wappinfo *w,
                          void *rawdata, float *spectra, float *lags,
                          fftwf_plan fftplan)
// This routine converts a single point of WAPP lags
// into a filterbank style spectrum array of floats.
// Van Vleck corrections are applied.
{
    int ii, ifnum = 0, wappnum = 0, index = 0;
    double power;

    // Loop over the WAPPs
    for (wappnum = 0; wappnum < numwapps; wappnum++) {
        index = wappnum * w->numifs * w->numchans;

        // Loop over the IFs
        for (ifnum = 0; ifnum < w->numifs; ifnum++, index += w->numchans) {

            // Fill lag array with scaled CFs
            if (w->bits_per_lag == 16) {
                unsigned short *sdata = (unsigned short *) rawdata;
                for (ii = 0; ii < w->numchans; ii++)
                    lags[ii] = w->corr_scale * sdata[ii + index] - 1.0;
            } else {
                unsigned int *idata = (unsigned int *) rawdata;
                for (ii = 0; ii < w->numchans; ii++)
                    lags[ii] = w->corr_scale * idata[ii + index] - 1.0;
            }

            // Calculate power
            power = inv_cerf(lags[0]);
            power = 0.1872721836 / (power * power);

            // Apply Van Vleck Corrections to the Lags
            if (w->corr_level == 3)
                vanvleck3lev(lags, w->numchans);
            else
                vanvleck9lev(lags, w->numchans);

            for (ii = 0; ii < w->numchans; ii++)
                lags[ii] *= power;

            // FFT the ACF lags (which are real and even) -> real and even FFT
            fftwf_execute(fftplan);

#if 0
            printf("\n");
            for (ii = 0; ii < w->numchans; ii++)
                printf("%d  %.7g\n", ii, lags[ii]);
            printf("\n");
            exit(0);
#endif

            // Reverse the band if it needs it
            if (w->invertband) {
                float tempzz = 0.0, *loptr, *hiptr;
                loptr = lags + 0;
                hiptr = lags + w->numchans - 1;
                for (ii = 0; ii < w->numchans / 2; ii++, loptr++, hiptr--) {
                    SWAP(*loptr, *hiptr);
                }
            }
            // Copy the spectra into the proper portion of the spectra array
            for (ii = 0; ii < w->numchans; ii++)
                spectra[index + ii] = lags[ii];
        }
    }
}



static double inv_cerf(double input)
/* Approximation for Inverse Complementary Error Function */
{
    static double numerator_const[3] = {
        1.591863138, -2.442326820, 0.37153461
    };
    static double denominator_const[3] = {
        1.467751692, -3.013136362, 1.0
    };
    double num, denom, temp_data, temp_data_srq, erf_data;

    erf_data = 1.0 - input;
    temp_data = erf_data * erf_data - 0.5625;
    temp_data_srq = temp_data * temp_data;
    num = erf_data * (numerator_const[0] +
                      (temp_data * numerator_const[1]) +
                      (temp_data_srq * numerator_const[2]));
    denom = denominator_const[0] + temp_data * denominator_const[1] +
        temp_data_srq * denominator_const[2];
    return num / denom;
}


#define NO    0
#define YES   1
/*------------------------------------------------------------------------*
 * Van Vleck Correction for 9-level sampling/correlation
 *  Samples {-4,-3,-2,-1,0,1,2,3,4}
 * Uses Zerolag to adjust correction
 *   data_array -> Points into ACF of at least 'count' points
 * This routine takes the first value as the zerolag and corrects the
 * remaining count-1 points.  Zerolag is set to a normalized 1
 * NOTE - The available routine works on lags normaized to -16<rho<16, so
 *  I need to adjust the values before/after the fit
 * Coefficent ranges
 *   c1
 *     all
 *   c2
 *     r0 > 4.5
 *     r0 < 2.1
 *     rest
 * NOTE - correction is done INPLACE ! Original values are destroyed
 * As reported by M. Lewis -> polynomial fits are OK, but could be improved
 *------------------------------------------------------------------------*/
static void vanvleck9lev(float *rho, int npts)
{
    double acoef[5], dtmp, zl;
    int i;
    static double coef1[5] =
        { 1.105842267, -0.053258115, 0.011830276, -0.000916417, 0.000033479 };
    static double coef2rg4p5[5] =
        { 0.111705575, -0.066425925, 0.014844439, -0.001369796, 0.000044119 };
    static double coef2rl2p1[5] =
        { 1.285303775, -1.472216011, 0.640885537, -0.123486209, 0.008817175 };
    static double coef2rother[5] =
        { 0.519701391, -0.451046837, 0.149153116, -0.021957940, 0.001212970 };
    static double coef3rg2p0[5] =
        { 1.244495105, -0.274900651, 0.022660239, -0.000760938, -1.993790548 };
    static double coef3rother[5] =
        { 1.249032787, 0.101951346, -0.126743165, 0.015221707, -2.625961708 };
    static double coef4rg3p15[5] =
        { 0.664003237, -0.403651682, 0.093057131, -0.008831547, 0.000291295 };
    static double coef4rother[5] =
        { 9.866677289, -12.858153787, 6.556692205, -1.519871179, 0.133591758 };
    static double coef5rg4p0[4] =
        { 0.033076469, -0.020621902, 0.001428681, 0.000033733 };
    static double coef5rg2p2[4] =
        { 5.284269565, 6.571535249, -2.897741312, 0.443156543 };
    static double coef5rother[4] =
        { -1.475903733, 1.158114934, -0.311659264, 0.028185170 };

    zl = rho[0] * 16;
    /* ro = *rho;                 */
    /*  for(i=0; i<npts; i++)     */
    /*    (rho+i) *= *(rho+i)/ro; */
    acoef[0] =
        ((((coef1[4] * zl + coef1[3]) * zl + coef1[2]) * zl +
          coef1[1]) * zl + coef1[0]);
    if (zl > 4.50)
        acoef[1] =
            ((((coef2rg4p5[4] * zl + coef2rg4p5[3]) * zl +
               coef2rg4p5[2]) * zl + coef2rg4p5[1]) * zl + coef2rg4p5[0]);
    else if (zl < 2.10)
        acoef[1] =
            ((((coef2rl2p1[4] * zl + coef2rl2p1[3]) * zl +
               coef2rl2p1[2]) * zl + coef2rl2p1[1]) * zl + coef2rl2p1[0]);
    else
        acoef[1] =
            ((((coef2rother[4] * zl + coef2rother[3]) * zl +
               coef2rother[2]) * zl + coef2rother[1]) * zl + coef2rother[0]);
    if (zl > 2.00)
        acoef[2] =
            coef3rg2p0[4] / zl +
            (((coef3rg2p0[3] * zl + coef3rg2p0[2]) * zl +
              coef3rg2p0[1]) * zl + coef3rg2p0[0]);
    else
        acoef[2] =
            coef3rother[4] / zl +
            (((coef3rother[3] * zl + coef3rother[2]) * zl +
              coef3rother[1]) * zl + coef3rother[0]);
    if (zl > 3.15)
        acoef[3] =
            ((((coef4rg3p15[4] * zl + coef4rg3p15[3]) * zl +
               coef4rg3p15[2]) * zl + coef4rg3p15[1]) * zl + coef4rg3p15[0]);
    else
        acoef[3] =
            ((((coef4rg3p15[4] * zl + coef4rother[3]) * zl +
               coef4rother[2]) * zl + coef4rother[1]) * zl + coef4rother[0]);
    if (zl > 4.00)
        acoef[4] =
            (((coef5rg4p0[3] * zl + coef5rg4p0[2]) * zl +
              coef5rg4p0[1]) * zl + coef5rg4p0[0]);
    else if (zl < 2.2)
        acoef[4] =
            (((coef5rg2p2[3] * zl + coef5rg2p2[2]) * zl +
              coef5rg2p2[1]) * zl + coef5rg2p2[0]);
    else
        acoef[4] =
            (((coef5rother[3] * zl + coef5rother[2]) * zl +
              coef5rother[1]) * zl + coef5rother[0]);
    for (i = 1; i < npts; i++) {
        dtmp = rho[i];
        rho[i] =
            ((((acoef[4] * dtmp + acoef[3]) * dtmp + acoef[2]) * dtmp +
              acoef[1]) * dtmp + acoef[0]) * dtmp;
    }
    rho[0] = 1.0;
    return;
}

/*------------------------------------------------------------------------*
 * Van Vleck Correction for 3-level sampling/correlation
 *  Samples {-1,0,1}
 * Uses Zerolag to adjust correction
 *   data_array -> Points into ACF of at least 'count' points
 * This routine takes the first value as the zerolag and corrects the
 * remaining count-1 points.  Zerolag is set to a normalized 1
 *
 * NOTE - correction is done INPLACE ! Original values are destroyed
 *------------------------------------------------------------------------*/
static void vanvleck3lev(float *rho, int npts)
{
    double lo_u[3], lo_h[3];
    double high_u[5], high_h[5];
    double lo_coefficient[3];
    double high_coefficient[5];
    double zho, zho_3;
    double temp_data;
    double temp_data_1;
    int ichan, ico, flag_any_high;
    static double lo_const[3][4] = {
        {0.939134371719, -0.567722496249, 1.02542540932, 0.130740914912},
        {-0.369374472755, -0.430065136734, -0.06309459132, -0.00253019992917},
        {0.888607422108, -0.230608118885, 0.0586846424223, 0.002012775510695}
    };
    static double high_const[5][4] = {
        {-1.83332160595, 0.719551585882, 1.214003774444, 7.15276068378e-5},
        {1.28629698818, -1.45854382672, -0.239102591283, -0.00555197725185},
        {-7.93388279993, 1.91497870485, 0.351469403030, 0.00224706453982},
        {8.04241371651, -1.51590759772, -0.18532022393, -0.00342644824947},
        {-13.076435520, 0.769752851477, 0.396594438775, 0.0164354218208}
    };

    /* Perform Lo correction on All data that is not flaged 
       for high correction  */
    zho = (double) rho[0];
    zho_3 = zho * zho * zho;
    lo_u[0] = zho;
    lo_u[1] = zho_3 - (61.0 / 512.0);
    lo_u[2] = zho - (63.0 / 128.0);
    lo_h[0] = zho * zho;
    lo_h[2] = zho_3 * zho_3 * zho;      /* zlag ^7 */
    lo_h[1] = zho * lo_h[2];    /* zlag ^8 */
    /* determine lo-correct coefficents - */
    for (ico = 0; ico < 3; ico++) {
        lo_coefficient[ico] =
            (lo_u[ico] *
             (lo_u[ico] *
              (lo_u[ico] * lo_const[ico][0] + lo_const[ico][1]) +
              lo_const[ico][2]) + lo_const[ico][3]) / lo_h[ico];
    }
    /* perform correction -- */
    for (ichan = 1, flag_any_high = NO; ichan < npts; ichan++) {
        temp_data = (double) rho[ichan];
        if (fabs(temp_data) > 0.199) {
            if (flag_any_high == NO) {
                high_u[0] = lo_h[2];    /* zlag ^7 */
                high_u[1] = zho - (63.0 / 128.0);
                high_u[2] = zho * zho - (31.0 / 128.0);
                high_u[3] = zho_3 - (61.0 / 512.0);
                high_u[4] = zho - (63.0 / 128.0);
                high_h[0] = lo_h[1];    /* zlag ^8 */
                high_h[1] = lo_h[1];    /* zlag ^8 */
                high_h[2] = lo_h[1] * zho_3 * zho;      /* zlag ^12 */
                high_h[3] = lo_h[1] * lo_h[1] * zho;    /* zlag ^17 */
                high_h[4] = high_h[3];  /* zlag ^17 */
                for (ico = 0; ico < 5; ico++) {
                    high_coefficient[ico] =
                        (high_u[ico] *
                         (high_u[ico] *
                          (high_u[ico] * high_const[ico][0] +
                           high_const[ico][1]) + high_const[ico][2]) +
                         high_const[ico][3]) / high_h[ico];
                }
                flag_any_high = YES;
            }
            temp_data_1 = fabs(temp_data * temp_data * temp_data);
            rho[ichan] =
                (temp_data *
                 (temp_data_1 *
                  (temp_data_1 *
                   (temp_data_1 *
                    (temp_data_1 * high_coefficient[4] +
                     high_coefficient[3]) + high_coefficient[2]) +
                   high_coefficient[1]) + high_coefficient[0]));
        } else {
            temp_data_1 = temp_data * temp_data;
            rho[ichan] =
                (temp_data *
                 (temp_data_1 *
                  (temp_data_1 * lo_coefficient[2] + lo_coefficient[1]) +
                  lo_coefficient[0]));
        }
    }
    rho[0] = 1.0;
}
