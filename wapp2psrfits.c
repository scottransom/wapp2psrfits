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
    float offset, scale, datum, packdatum;
    float *lags, *fspects, *datachunk;
    long long N = 0, totsize;
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
    printf("Code has git hash:  '%s'\n\n", GITHASH);

    if (!(cmd->numbits == 4 || cmd->numbits == 8)) {
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
    N = get_WAPP_info(cmd->argv[0], infiles, numfiles, cmd->numwapps, &hdr, &w);
    // Print a summary...
    print_WAPP_hdr(hdr);
    printf("\nFound a total of %lld samples.\n", N);
    // Are we explicitly flipping the band?
    if (cmd->invertP) {
        printf("\nForcing a flip of the band since '-i' is specified.\n");
        w.invertband = (w.invertband == 1) ? 0 : 1;
    }

    // Write a copy of the first(!) original WAPP header
    {
        long posn = ftell(infiles[0]);
        unsigned char *hdr = gen_bvect(w.header_size);
        char hdrname[200];
        FILE *hdrfile;

        // Go to the beginning of the first file
        chkfileseek(infiles[0], 0, 1, SEEK_SET);
        // Read the ASCII and binary headers together
        chkfread(hdr, w.header_size, 1, infiles[0]);
        // The following should be redundant as we should already be there
        // after reading the header, but just in case...
        chkfileseek(infiles[0], posn, 1, SEEK_SET);
        // Determine the output WAPP header file name
        sprintf(hdrname, "%s.wapp_hdr", cmd->outfile);
        printf("Copying the header from the first WAPP file to '%s'\n", hdrname);
        // Open, write the header, and close the header file
        hdrfile = chkfopen(hdrname, "wb");
        chkfwrite(hdr, w.header_size, 1, hdrfile);
        fclose(hdrfile);
        free(hdr);
    }

    // Prep the psrfits structure
    fill_psrfits_struct(cmd->numwapps, cmd->numbits, hdr, &w, &pf);
    close_parse(hdr);
    spec_per_row = pf.hdr.nsblk;
    numrows = N / spec_per_row;
    totsize = numrows * pf.sub.bytes_per_subint;
    if (totsize / 1073741824.0 < cmd->outlenGB) {
        pf.multifile = 0;
        pf.rows_per_file = numrows;
    } else {
        pf.multifile = 1;
        pf.rows_per_file = (int) ((cmd->outlenGB * 1073741824.0)
                                  / (double) pf.sub.bytes_per_subint);
    }
    printf("PSRFITS file(s) will have %d samples per row and %d rows.\n\n",
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
    // See:  http://www.cv.nrao.edu/~pdemores/wapp/
    lags = (float *) fftwf_malloc((w.numchans) * sizeof(float));
    fftplan = fftwf_plan_r2r_1d(w.numchans, lags,
                                lags, FFTW_REDFT01, FFTW_PATIENT);

    // Create the PSRFITS file
    strcpy(pf.basefilename, cmd->outfile);
    psrfits_create(&pf);

    // Loop over the data

    // The rows in the output file
    for (rownum = 0; rownum < numrows; rownum++) {
        printf("\rWorking on row %d", rownum + 1);
        fflush(stdout);

        // Loop over all the spectra per row
        for (specnum = 0; specnum < spec_per_row; specnum++) {
            read_WAPP_lags(infiles, numfiles, cmd->numwapps, rawdata, &w);
            WAPP_lags_to_spectra(cmd->numwapps, &w, rawdata,
                                 fspects + specnum * pf.hdr.nchan * pf.hdr.npol,
                                 lags, fftplan);
        }

        // Loop over all the channels:
        for (ichan = 0; ichan < pf.hdr.nchan * pf.hdr.npol; ichan++) {

            // Populate datachunk[] by picking out all time samples for ichan
            for (itsamp = 0; itsamp < spec_per_row; itsamp++) {
                datachunk[itsamp] =
                    fspects[ichan + itsamp * pf.hdr.nchan * pf.hdr.npol];
            }

            // Compute the statistics here, and put the offsets and scales in
            // pf.sub.dat_offsets[] and pf.sub.dat_scales[]

            if (rescale(datachunk, spec_per_row, cmd->numbits, &offset, &scale) != 0) {
                printf("Rescale routine failed!\n");
                return (-1);
            }
            pf.sub.dat_offsets[ichan] = offset;
            pf.sub.dat_scales[ichan] = scale;

            // Since we have the offset and scale ready, rescale the data:
            for (itsamp = 0; itsamp < spec_per_row; itsamp++) {
                datum = (scale == 0.0) ? 0.0 :
                    roundf((datachunk[itsamp] - offset) / scale);
                // Check that it lies between 0 and (2^Nbits -1)
                datum = (datum < 0) ? 0.0 : datum;
                if (cmd->numbits == 4) {
                    datum = (datum > 15.0) ? 15.0 : datum;
                } else if (cmd->numbits == 8) {
                    datum = (datum > 255.0) ? 255.0 : datum;
                } else {
                    printf("This can't be happening!\n");
                }
                fspects[ichan + itsamp * pf.hdr.nchan * pf.hdr.npol] = datum;
            }
            // Now fspects[ichan] contains rescaled and clipped floats.
        }

        // Then do the conversion to 4-bits or 8-bits and store the
        // results in pf.sub.data[] 
        if (cmd->numbits == 4) {
            for (itsamp = 0; itsamp < spec_per_row; itsamp++) {
                datidx = itsamp * pf.hdr.nchan * pf.hdr.npol;
                for (ichan = 0; ichan < pf.hdr.nchan * pf.hdr.npol;
                     ichan += 2, datidx += 2) {
                    packdatum = fspects[datidx] * 16 + fspects[datidx + 1];
                    pf.sub.data[datidx / 2] = (unsigned char) packdatum;
                }
            }
        } else {                // cmd->numbits==8
            for (itsamp = 0; itsamp < spec_per_row; itsamp++) {
                datidx = itsamp * pf.hdr.nchan * pf.hdr.npol;
                for (ichan = 0; ichan < pf.hdr.nchan * pf.hdr.npol;
                     ichan++, datidx++) {
                    //if (fspects[datidx] > 256.0 || fspects[datidx] < 0.0) {
                    //    printf("Yikes!  %d  %d  %.7g\n", itsamp, ichan, fspects[datidx]);
                    //}
                    pf.sub.data[datidx] = (unsigned char) fspects[datidx];
                }
            }
        }

        // Now write the row...
        pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;
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
