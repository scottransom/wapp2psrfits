#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fftw3.h"
#include "chkio.h"
#include "wapp.h"
#include "vectors.h"
#include "wapp2psrfits_cmd.h"
#include "psrfits.h"
#include "rescale.h"

extern short transpose_float(float *a, int nx, int ny, unsigned char *move,
                             int move_size);

int main(int argc, char *argv[])
{
  int numfiles, ii, numrows, rownum, ichan, itsamp, datidx;
    int spec_per_row, specnum, status, move_size;
    short trtn;
    float offset, scale, datum, packdatum; 
    float *lags, *fspects, *datachunk;
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
#ifdef DEBUG
    showOptionValues();
#endif
    
    printf("\n           WAPP to PSRFITs Conversion Code\n");
    printf("            by S. Ransom & S. Chatterjee\n\n");
    
    if (!(cmd->numbits==4 || cmd->numbits==8)) {
        printf("ERROR:  -b (# of output bits) must be 4 or 8!\n");
        exit(1);
    }

    numfiles = cmd->argc;
    infiles = (FILE **) malloc(numfiles * sizeof(FILE *));
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
    fill_psrfits_struct(cmd->numwapps, cmd->numbits, hdr, &w, &pf);
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
    // Staging array for the rescaling process: 
    datachunk = gen_fvect(spec_per_row);
    // (size is usually 1 sec / 64 microsec = 15625)

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
        for (specnum = 0 ; specnum < spec_per_row ; specnum++) {
            read_WAPP_lags(infiles, numfiles, cmd->numwapps, rawdata, &w);
            WAPP_lags_to_spectra(cmd->numwapps, &w, rawdata, 
                                 fspects + specnum * pf.hdr.nchan * pf.hdr.npol, 
                                 lags, fftplan);
        }

        // Loop over all the channels:
        for (ichan = 0 ; ichan < pf.hdr.nchan * pf.hdr.npol ; ichan++){

            // Populate datachunk[] by picking out all time samples for ichan
            for (itsamp = 0 ; itsamp < spec_per_row ; itsamp++){
                datachunk[itsamp] = 
                    fspects[ichan + itsamp * pf.hdr.nchan * pf.hdr.npol];
            }
            
            // Compute the statistics here, and put the offsets and scales in
            // pf.sub.dat_offsets[] and pf.sub.dat_scales[]
            
            if (rescale(datachunk, spec_per_row, cmd->numbits, &offset, &scale)!=0){
                printf("Rescale routine failed!\n");
                return(-1);
            }
            pf.sub.dat_offsets[ichan] = offset;
            pf.sub.dat_scales[ichan] = scale;
            
            // Since we have the offset and scale ready, rescale the data:
            for (itsamp = 0 ; itsamp < spec_per_row ; itsamp++){
                datum = (scale==0.0) ? 0.0 : \
                    roundf((datachunk[itsamp] - offset) / scale);
                fspects[ichan + itsamp * pf.hdr.nchan * pf.hdr.npol] = datum;
            }	  
            // Now fspects[ichan] contains rescaled floats.
        }

        // Then do the conversion to 4-bits or 8-bits and store the
        // results in pf.sub.data[] 
        if (cmd->numbits==4) {
            for (itsamp = 0 ; itsamp < spec_per_row ; itsamp++){
                datidx = itsamp * pf.hdr.nchan * pf.hdr.npol;
                for (ichan=0 ; ichan < pf.hdr.nchan * pf.hdr.npol ; 
                     ichan+=2, dataidx++){
                    packdatum = fspects[datidx] * 16 + fspects[datidx + 1];
                    pf.sub.data[datidx/2] = (unsigned char)packdatum;
                }
            }
        } else {  // cmd->numbits==8
            for (itsamp = 0 ; itsamp < spec_per_row ; itsamp++){
                datidx = itsamp * pf.hdr.nchan * pf.hdr.npol;
                for (ichan=0 ; ichan < pf.hdr.nchan * pf.hdr.npol ; 
                     ichan++, dataidx++){
                    //if (fspects[datidx] > 256.0 || fspects[datidx] < 0.0) {
                    //    printf("Yikes!  %d  %d  %.7g\n", itsamp, ichan, fspects[datidx]);
                    //}
                    pf.sub.data[datidx] = (unsigned char)fspects[datidx];
                }
            }
        }
	    
        // Now write the row...
        pf.sub.offs = (pf.tot_rows+0.5) * pf.sub.tsubint;
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
    free(datachunk);
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
