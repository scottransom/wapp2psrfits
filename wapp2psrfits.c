#include <stdio.h>
#include <stdlib.h>
#include "chkio.h"
#include "wapp.h"
#include "wapp2psrfits_cmd.h"
#include "psrfits.h"

int main(int argc, char *argv[])
{
    int numfiles, ii;
    long long N=0;
    FILE **infiles;
    struct HEADERP hdr;
    struct wappinfo w;
    Cmdline *cmd;
    
    if (argc == 1) {
        Program = argv[0];
        usage();
        exit(1);
    }

    /* Parse the command line using the excellent program Clig */
    
    cmd = parseCmdline(argc, argv);
    numfiles = cmd->argc;
    infiles = (FILE **) malloc(numfiles * sizeof(FILE *));
    
#ifdef DEBUG
    showOptionValues();
#endif
    
    printf("\n           WAPP to PSRFITs Conversion Code\n");
    printf("              by S. Ransom & S. Chatterjee\n\n");
    
    // Open the input files
    printf("Reading input data from:\n");
    for (ii = 0; ii < numfiles; ii++) {
        printf("  '%s'\n", cmd->argv[ii]);
        infiles[ii] = chkfopen(cmd->argv[ii], "rb");
    }

    // Get the basic information from them
    N = get_WAPP_info(infiles, numfiles, cmd->numwapps, &hdr, &w);
    printf("Found a total of %lld samples.\n", N);
                      
        
    /* The following should be freed sometime... */
    //lags = (float *) fftwf_malloc((w->numchans + 1) * sizeof(float));
    //fftplan = fftwf_plan_r2r_1d(w->numchans + 1, lags, 
    //                            lags, FFTW_REDFT00, FFTW_PATIENT);
    
    return 0;
}
