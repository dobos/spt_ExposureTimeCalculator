#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "libgsetc.h"

/* Compute the line spread function (LSF) */
/* Written by L.Dobos base on the original ETC */

typedef struct {
    char* spectroConfig;
    char* psfFile;
    int N;
} PARAMS;

PARAMS getParams(int argc, char* argv[]) {
    PARAMS params = {
        .spectroConfig = gsGetArgPositional(argc, argv, 1),
        .psfFile = gsGetArgPositional(argc, argv, 2),
        .N = atoi(gsGetArgPositional(argc, argv, 3))
    };
    return params;
}

int main(int argc, char* argv[]) {
    SPECTRO_ATTRIB spectro;

    int i_arm, i;
    long ip;
    double lambda;
    double tr;
    double* fr;
    FILE* fp;

    gsPrintCompilerFlags();

    /* Process command-line arguments */
    PARAMS params = getParams(argc, argv);

    /* Load config, do not degrade detector throughput */
    fp = gsOpenConfigFile(params.spectroConfig);
    gsReadSpectrographConfig(fp, &spectro, 1.0);
    gsCloseConfigFile(fp);

    /* Pretabulate the point spread function
     * It depends on the wavelength and the distance (in pixels) from the center 
     * Compute it in 2N + 1 pixels, centered on every pixel */

    fp = gsOpenOutputFile(params.psfFile);
    
    fr = malloc((2 * params.N + 1) * sizeof(double));
    for (i_arm = 0; i_arm < spectro.N_arms; i_arm++) {
        for (ip = 0; ip < spectro.npix[i_arm]; ip++) {
            // Take wavelength at the center of the pixel
            lambda = spectro.lmin[i_arm] + (ip + 0.5) * spectro.dl[i_arm];   
            gsSpectroDist(&spectro, i_arm, lambda, params.N, 0, 2 * params.N + 1, fr, 1);
            tr = gsFracTrace(&spectro, i_arm, lambda, 0);
            fprintf(fp, "%d %ld %f ", spectro_arm(&spectro, i_arm), ip, lambda);
            for (i = 0; i < 2 * params.N + 1; i ++) {
                fprintf(fp, "%f ", fr[i]);
            }
            fprintf(fp, "%f", tr);
            fprintf(fp, "\n");
        }
    }
    free(fr);

    gsCloseOutputFile(fp);
}