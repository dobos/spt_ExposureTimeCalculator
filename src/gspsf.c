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
    int oversampling;
    int binary;
} PARAMS;

PARAMS getParams(int argc, char* argv[]) {
    PARAMS params = {
        .spectroConfig = gsGetArgPositional(argc, argv, 1),
        .psfFile = gsGetArgPositional(argc, argv, 2),
        .N = atoi(gsGetArgPositional(argc, argv, 3)),
        .oversampling = gsGetArgNamedInt(argc, argv, "--oversampling", 1),
        .binary = gsGetArgNamedBoolean(argc, argv, "--binary")
    };
    return params;
}

int main(int argc, char* argv[]) {
    SPECTRO_ATTRIB spectro;

    int i_arm, i, arm;
    long ip;
    double lambda;
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
    
    fr = malloc(params.N * params.oversampling * sizeof(double));
    for (i_arm = 0; i_arm < spectro.N_arms; i_arm++) {
        for (ip = 0; ip < spectro.npix[i_arm]; ip++) {
            arm = spectro_arm(&spectro, i_arm);
            // Take wavelength at the center of the pixel
            lambda = spectro.lmin[i_arm] + (ip + 0.5) * spectro.dl[i_arm];   
            gsSpectroDist(&spectro, i_arm, lambda, 0.5 * params.N, 0, params.N, fr, params.oversampling);

            if (params.binary) {
                fwrite(fr, sizeof(double), params.N * params.oversampling, fp);
            } else {
                fprintf(fp, "%d %ld %f ", arm, ip, lambda);
                for (i = 0; i < params.N * params.oversampling; i ++) {
                    fprintf(fp, "%f ", fr[i]);
                }
                fprintf(fp, "\n");
            }
        }
    }
    free(fr);

    gsCloseOutputFile(fp);
}