#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fftw3.h"
#include "chkio.h"
#include "wapp.h"
#include "vectors.h"
#include "wapp2psrfits_cmd.h"
#include "psrfits.h"

extern short transpose_float(float *a, int nx, int ny, unsigned char *move,
                             int move_size);

int main(int argc, char *argv[])
{
    int numfiles, ii, numrows, rownum;
    int spec_per_row, specnum, status, move_size;
    short trtn;
    float *lags, *fspects;
    long long N=0;
    unsigned char *move, *rawdata;
    FILE **infiles;
    struct HEADERP *hdr;
    struct wappinfo w;
    struct psrfits pf;
    Cmdline *cmd;
    fftwf_plan fftplan;
    
    if (argc == 1) {
        Program = argv[0];
        usage();
        exit(1);
    }

    // Parse the command line using the excellent program Clig
    cmd = parseCmdline(argc, argv);
    numfiles = cmd->argc;
    infiles = (FILE **) malloc(numfiles * sizeof(FILE *));
    
#ifdef DEBUG
    showOptionValues();
#endif
    
    printf("\n           WAPP to PSRFITs Conversion Code\n");
    printf("            by S. Ransom & S. Chatterjee\n\n");
    
    // Open the input files
    printf("Reading input data from:\n");
    for (ii = 0; ii < numfiles; ii++) {
        printf("  '%s'\n", cmd->argv[ii]);
        infiles[ii] = chkfopen(cmd->argv[ii], "rb");
    }

    // Get the basic information from them
    N = get_WAPP_info(infiles, numfiles, cmd->numwapps, &hdr, &w);
    printf("Found a total of %lld samples.\n", N);
                      
    // Prep the psrfits structure
    fill_psrfits_struct(cmd->numwapps, hdr, &w, &pf);
    close_parse(hdr);
    spec_per_row = pf.hdr.nsblk;
    numrows = N / spec_per_row;
    pf.rows_per_file = numrows; // NOTE: this will make a _single_ file!  FIXME!
    printf("PSRFITS will have %d samples per row and %d rows.\n", 
           spec_per_row, numrows);

    // Create the arrays we will need
    //
    // These are the raw lags for a single sample from all WAPPs
    rawdata = gen_bvect(w.bytes_per_sample * cmd->numwapps);
    // These are the floating-point converted spectra for a full row
    fspects = gen_fvect(cmd->numwapps * w.numifs * w.numchans * spec_per_row);
    // For the transposes
    move_size = (spec_per_row + pf.hdr.nchan * pf.hdr.npol) / 2;
    move = gen_bvect(move_size);

    // Create the FFTW plan
    lags = (float *) fftwf_malloc((w.numchans + 1) * sizeof(float));
    fftplan = fftwf_plan_r2r_1d(w.numchans + 1, lags, 
                                lags, FFTW_REDFT00, FFTW_PATIENT);

    // Create the PSRFITS file
    strcpy(pf.basefilename, cmd->outfile);
    psrfits_create(&pf);

    // Loop over the data

    // The rows in the output file
    for (rownum = 0 ; rownum < numrows ; rownum++) {
        printf("\rWorking on row %d", rownum+1);
        fflush(stdout);

        // Loop over all the spectra per row
        for (specnum = 0 ; spec_per_row < numrows ; specnum++) {
            read_WAPP_lags(infiles, numfiles, cmd->numwapps, rawdata, &w);
            WAPP_lags_to_spectra(cmd->numwapps, &w, rawdata, 
                                 fspects + specnum * pf.hdr.nchan * pf.hdr.npol, 
                                 lags, fftplan);
        }

        // Transpose the array so that it is grouped by channels 
        // rather than spectra
        if ((trtn = transpose_float(fspects, spec_per_row, 
                                    pf.hdr.nchan * pf.hdr.npol,
                                    move, move_size)) < 0) {
            printf("\nError (%d) in transpose_float().\n\n", trtn);
            exit(1);
        }
        
        // Compute the statistics here, and put the offsets and scales in
        // pf.sub.dat_offsets[] and pf.sub.dat_scales[]

        // Then do the conversion to 4- or 8-bits and store the
        // results in pf.sub.data[], transposing as we go...

        // Now write the row...
        status = psrfits_write_subint(&pf);
        if (status) {
            printf("\nError (%d) writing PSRFITS...\n\n", status);
            break;
        }
    }
    printf("\n");

    // Close the PSRFITS file
    psrfits_close(&pf);

    // Free the structure arrays too...
    free(pf.sub.dat_freqs);
    free(pf.sub.dat_weights);
    free(pf.sub.dat_offsets);
    free(pf.sub.dat_scales);
    free(pf.sub.data);
    free(fspects);
    free(rawdata);
    free(lags);
    free(move);
    return 0;
}
