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
    // int ia;
    // for (ia = 0; ia < spectro->N_arms; ia ++) {
    //    spectro->npix[ia] = 128;
    //    spectro->lmax[ia] = spectro->lmin[ia] * 128 * spectro->dl[ia];
    //}

    spectro->diffuse_stray=0.2;                     // Bump this up a bit to see if works
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

void ASSERT(int expr, int line) {
    if (!expr) {
        fprintf(stderr, "Assertion failed at line %d", line);
        exit(1);
    }
}

/* Pretabulate the point spread function at the center of each pixel
 * */
void test_SpectroDist(FILE* fp, SPECTRO_ATTRIB *spectro, int N) {
    printf("Running test %s.\n", __func__);

    int i_arm, i;
    long ip;
    double lambda;
    double* fr;
    
    fr = malloc(N * sizeof(double));
    for (i_arm = 0; i_arm < spectro->N_arms; i_arm++) {
        for (ip = 0; ip < spectro->npix[i_arm]; ip++) {
            lambda = gsGetLambdaCenter(spectro, i_arm, ip);
            gsSpectroDist(spectro, i_arm, lambda, N / 2.0, 0, N, fr, 1);
            fprintf(fp, "%d %ld %f ", spectro_arm(spectro, i_arm), ip, lambda);
            for (i = 0; i < N; i ++) {
                fprintf(fp, "%f ", fr[i]);
            }
            fprintf(fp, "\n");

            for (i = 0; i < N / 2; i ++) {
                ASSERT(fabs(fr[i] - fr[N - i - 1]) < 1e-7, __LINE__);
            }
        }
    }
    free(fr);
}

/* Pretabulate the point spread function at an offset from the pixel center
 * */
void test_SpectroDist_Oversampled(FILE* fp, SPECTRO_ATTRIB *spectro, int N) {
    printf("Running test %s.\n", __func__);

    int i_arm, i;
    long ip;
    double lambda;
    double* fr;
    
    fr = malloc(N * spectro->oversampling * sizeof(double));
    for (i_arm = 0; i_arm < spectro->N_arms; i_arm++) {
        for (ip = 0; ip < spectro->npix[i_arm]; ip++) {
            lambda = gsGetLambdaCenter(spectro, i_arm, ip);
            gsSpectroDist(spectro, i_arm, lambda, N / 2.0, 0, N, fr, spectro->oversampling);
            fprintf(fp, "%d %ld %f ", spectro_arm(spectro, i_arm), ip, lambda);
            for (i = 0; i < N * spectro->oversampling; i ++) {
                fprintf(fp, "%f ", fr[i]);
            }
            fprintf(fp, "\n");
        }
    }
    free(fr);
}

void test_FracTrace(FILE* fp, SPECTRO_ATTRIB *spectro) {
    printf("Running test %s.\n", __func__);

    int i_arm;
    long ip;
    double lambda, frac_0, frac_1;

    for (i_arm = 0; i_arm < spectro->N_arms; i_arm++) {
        for (ip = 0; ip < spectro->npix[i_arm]; ip++) {
            lambda = gsGetLambdaCenter(spectro, i_arm, ip);
            frac_0 = gsFracTrace(spectro, i_arm, lambda, 0);
            frac_1 = gsFracTrace(spectro, i_arm, lambda, 1);
            fprintf(fp, "%d %ld %f ", spectro_arm(spectro, i_arm), ip, lambda);
            fprintf(fp, "%f %f\n", frac_0, frac_1);
        }
    }
}

void test_AtmTrans(FILE* fp, SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs) {
    printf("Running test %s.\n", __func__);

    int i_arm;
    long ip;
    double lambda, atm_0, atm_1, atm_2;

    for (i_arm = 0; i_arm < spectro->N_arms; i_arm++) {
        for (ip = 0; ip < spectro->npix[i_arm]; ip++) {
            lambda = gsGetLambdaCenter(spectro, i_arm, ip);
            atm_0 = gsAtmTrans(obs, lambda);
            atm_1 = gsAtmTransInst(spectro, obs, i_arm, lambda);
            atm_2 = gsAtmContOp(obs, lambda);
            fprintf(fp, "%d %ld %f ", spectro_arm(spectro, i_arm), ip, lambda);
            fprintf(fp, "%f %f %f\n", atm_0, atm_1, atm_2);
        }
    }
}

void test_GetSignal(FILE* fp, SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs) {
    printf("Running test %s.\n", __func__);

    int i_arm;
    long ip;
    double lambda;
    double *signal_1, *signal_2;

    for (i_arm = 0; i_arm < spectro->N_arms; i_arm++) {
        signal_1 = malloc(spectro->npix[i_arm] * sizeof(double));
        signal_2 = malloc(spectro->npix[i_arm] * sizeof(double));
        gsGetSignal(spectro, obs, i_arm, 680, 1e-17, 30, signal_1);
        gsGetSignal(spectro, obs, i_arm, spectro->lmin[i_arm] + 5 * spectro->dl[i_arm], 1e-17, 30, signal_2);
        for (ip = 0; ip < spectro->npix[i_arm]; ip++) {
            lambda = gsGetLambdaCenter(spectro, i_arm, ip);
            fprintf(fp, "%d %ld %f ", spectro_arm(spectro, i_arm), ip, lambda);
            fprintf(fp, "%f %f\n", signal_1[ip], signal_2[ip]);
        }
        free(signal_1);
        free(signal_2);
    }
}

void test_AddSkyLines(FILE* fp, SPECTRO_ATTRIB* spectro, OBS_ATTRIB* obs) {
    printf("Running test %s.\n", __func__);

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
    printf("Running test %s.\n", __func__);

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

void test_AddStrayLight(FILE* fp, SPECTRO_ATTRIB* spectro, OBS_ATTRIB* obs) {
    printf("Running test %s.\n", __func__);

    int ip;
    int i_arm ;
    double lambda, sample_factor;
    double **counts, **sky;

    gsAllocArmVectors(spectro, &counts);
    gsAllocArmVectors(spectro, &sky);
    for (i_arm = 0; i_arm < spectro->N_arms; i_arm++) {
        sample_factor = gsGetSampleFactor(spectro, i_arm);
        gsAddSkyContinuum(spectro, obs, i_arm, counts[i_arm]);
        gsAddLunarContinuum(spectro, obs, i_arm, counts[i_arm]);
        for (ip = 0; ip < spectro->npix[i_arm] * spectro->oversampling; ip++) {
            sky[i_arm][ip] = counts[i_arm][ip] / sample_factor;
            counts[i_arm][ip] = 0;
        }
        gsAddStrayLight(spectro, obs, i_arm, counts[i_arm], sky[i_arm]);
        for (ip = 0; ip < spectro->npix[i_arm] * spectro->oversampling; ip++) {
            lambda = spectro->lmin[i_arm] + (ip + 0.5) * spectro->dl[i_arm] / spectro->oversampling;
            fprintf(fp, "%d %d %f %f\n", spectro_arm(spectro, i_arm), ip, lambda, counts[i_arm][ip]);
        }
    }
    gsFreeArmVectors(spectro, counts);
    gsFreeArmVectors(spectro, sky);
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

    fp = gsOpenOutputFile("test_SpectroDist_Even.dat");
    test_SpectroDist(fp, &spectro, 4);
    gsCloseOutputFile(fp);

    fp = gsOpenOutputFile("test_SpectroDist_Odd.dat");
    test_SpectroDist(fp, &spectro, 5);
    gsCloseOutputFile(fp);

    fp = gsOpenOutputFile("test_SpectroDist_Oversampled_Even.dat");
    test_SpectroDist_Oversampled(fp, &spectro, 2);
    gsCloseOutputFile(fp);

    fp = gsOpenOutputFile("test_SpectroDist_Oversampled_Odd.dat");
    test_SpectroDist_Oversampled(fp, &spectro, 3);
    gsCloseOutputFile(fp);

    fp = gsOpenOutputFile("test_FracTrace.dat");
    test_FracTrace(fp, &spectro);
    gsCloseOutputFile(fp);

    fp = gsOpenOutputFile("test_AtmTrans.dat");
    test_AtmTrans(fp, &spectro, &obs);
    gsCloseOutputFile(fp);

    fp = gsOpenOutputFile("test_GetSignal.dat");
    test_GetSignal(fp, &spectro, &obs);
    gsCloseOutputFile(fp);

    fp = gsOpenOutputFile("test_AddSkyLines.dat");
    test_AddSkyLines(fp, &spectro, &obs);
    gsCloseOutputFile(fp);

    fp = gsOpenOutputFile("test_AddSkyContinuum.dat");
    test_AddSkyContinuum(fp, &spectro, &obs);
    gsCloseOutputFile(fp);

    fp = gsOpenOutputFile("test_AddStrayLight.dat");
    test_AddStrayLight(fp, &spectro, &obs);
    gsCloseOutputFile(fp);
}