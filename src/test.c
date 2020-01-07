#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "libgsetc.h"

typedef struct {
    char* spectroConfig;
    int oversampling;
} PARAMS;

PARAMS getParams(int argc, char* argv[]) {
    PARAMS params = {
        .spectroConfig = gsGetArgPositional(argc, argv, 1),
        .oversampling = atoi(gsGetArgPositional(argc, argv, 2))
    };
    return params;
}

void updateSpectroConfig(SPECTRO_ATTRIB* spectro, PARAMS* params) {
    /* Reduce number of pixels for faster testing */
    int ia;
    for (ia = 0; ia < spectro->N_arms; ia ++) {
        spectro->npix[ia] = 128;
        spectro->lmax[ia] = spectro->lmin[ia] * 128 * spectro->dl[ia];
    }

    spectro->oversampling = params->oversampling;
}

OBS_ATTRIB createObsConfig() {
    OBS_ATTRIB obs = {
        .seeing_fwhm_800 = 0.8,
        .elevation = 4000,
        .zenithangle = 0.0,
        .lunarphase = 0.0,
        .lunarangle = 180,
        .lunarZA = 180,
        .EBV = 0,
        .skytype = 0x00011005,
        .fieldangle = 0.0,
        .decent = 0.0,
        .t_exp = 450,
        .n_exp = 1,
        .r_eff = 0,
        .flags = 0
    };
    return obs;
}

/* Pretabulate the point spread function
 * It depends on the wavelength and the distance (in pixels) from the center
 * */
void test_LineSpreadFunction(FILE* fp, SPECTRO_ATTRIB *spectro, int N) {
    printf("Running test test_LineSpreadFunction");

    int i_arm, i;
    long ip;
    double lambda;
    double* fr;
    
    fr = malloc(sizeof(double) * N);
    for (i_arm = 0; i_arm < spectro->N_arms; i_arm++) {
        for (ip = 0; ip < spectro->npix[i_arm]; ip++) {
            lambda = spectro->lmin[i_arm] + (ip + 0.5) * spectro->dl[i_arm];
            gsSpectroDist(spectro, i_arm, lambda, 0, 0, N, fr);
            fprintf(fp, "%d %ld %f ", spectro_arm(spectro, i_arm), ip, lambda);
            for (i = 0; i < N; i ++) {
                fprintf(fp, "%f ", fr[i]);
            }
            fprintf(fp, "\n");
        }
    }
    free(fr);
}

void test_LineSpreadFunction_Oversampled(FILE* fp, SPECTRO_ATTRIB *spectro, int N) {
    printf("Running test test_LineSpreadFunction_Oversampled");

    int i_arm, i;
    long ip;
    double lambda;
    double* fr;
    
    fr = malloc(sizeof(double) * N);
    for (i_arm = 0; i_arm < spectro->N_arms; i_arm++) {
        for (ip = 0; ip < spectro->npix[i_arm]; ip++) {
            lambda = spectro->lmin[i_arm] + (ip + 0.5) * spectro->dl[i_arm];   
            gsSpectroDist_Oversampled(spectro, i_arm, lambda, 0, 0, N, fr);
            fprintf(fp, "%d %ld %f ", spectro_arm(spectro, i_arm), ip, lambda);
            for (i = 0; i < N; i ++) {
                fprintf(fp, "%f ", fr[i]);
            }
            fprintf(fp, "\n");
        }
    }
    free(fr);
}

void test_AddSkyLines(FILE* fp, SPECTRO_ATTRIB* spectro, OBS_ATTRIB* obs) {
    int ip;
    int i_arm ;
    double lambda;
    double** counts;

    gsAllocArmVectors(spectro, &counts);
    for (i_arm = 0; i_arm < spectro->N_arms; i_arm++) {
        gsAddSkyLines(spectro, obs, i_arm, counts[i_arm]);
        for (ip = 0; ip < spectro->npix[i_arm] * spectro->oversampling; ip++) {
            lambda = spectro->lmin[i_arm] + (ip + 0.5) * spectro->dl[i_arm] / spectro->oversampling;
            fprintf(fp, "%d %d %f %f\n", spectro_arm(spectro, i_arm), ip, lambda, counts[i_arm][ip]);
        }
    }
    gsFreeArmVectors(spectro, counts);
}

void test_AddSkyContinuum(FILE* fp, SPECTRO_ATTRIB* spectro, OBS_ATTRIB* obs) {
    int ip;
    int i_arm ;
    double lambda;
    double** counts;

    gsAllocArmVectors(spectro, &counts);
    for (i_arm = 0; i_arm < spectro->N_arms; i_arm++) {
        gsAddSkyContinuum(spectro, obs, i_arm, counts[i_arm]);
        gsAddLunarContinuum(spectro, obs, i_arm, counts[i_arm]);
        for (ip = 0; ip < spectro->npix[i_arm] * spectro->oversampling; ip++) {
            lambda = spectro->lmin[i_arm] + (ip + 0.5) * spectro->dl[i_arm] / spectro->oversampling;
            fprintf(fp, "%d %d %f %f\n", spectro_arm(spectro, i_arm), ip, lambda, counts[i_arm][ip]);
        }
    }
    gsFreeArmVectors(spectro, counts);
}

int main(int argc, char* argv[]) {
    FILE* fp;

    gsPrintCompilerFlags();

    /* Process command-line arguments */
    PARAMS params = getParams(argc, argv);

    /* Load config, do not degrade detector throughput */
    SPECTRO_ATTRIB spectro;
    fp = gsOpenConfigFile(params.spectroConfig);
    gsReadSpectrographConfig(fp, &spectro, 1.0);
    gsCloseConfigFile(fp);
    updateSpectroConfig(&spectro, &params);

    OBS_ATTRIB obs = createObsConfig();

    /* Run tests */

    /*
    fp = gsOpenOutputFile("test_LineSpreadFunction.dat");
    test_LineSpreadFunction(fp, &spectro, 8);
    gsCloseOutputFile(fp);
    */

    /*
    fp = gsOpenOutputFile("test_LineSpreadFunction_Oversampled.dat");
    test_LineSpreadFunction_Oversampled(stdout, &spectro, 8);
    gsCloseOutputFile(fp);
    */

    /*
    fp = gsOpenOutputFile("test_AddSkyLines.dat");
    test_AddSkyLines(fp, &spectro, &obs);
    gsCloseOutputFile(fp);
    */

    fp = gsOpenOutputFile("test_AddSkyContinuum.dat");
    test_AddSkyContinuum(fp, &spectro, &obs);
    gsCloseOutputFile(fp);
}