#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "libgsetc.h"

/* Compute sky model and systematic error */
/* Written by L.Dobos base on the original ETC */

typedef struct {
    char* spectroConfig;
    char* obsConfig;
    char* noiseFile;
    int addSkyLines;
    int addSkyContinuum;
    int addLunarContinuum;
    int addStrayLight;
    int addDarkNoise;
    int addReadoutNoise;
} PARAMS;

PARAMS getParams(int argc, char* argv[]) {
    PARAMS params = {
        .spectroConfig = gsGetArgPositional(argc, argv, 1),
        .obsConfig = gsGetArgPositional(argc, argv, 2),
        .noiseFile = gsGetArgPositional(argc, argv, 3),
        .addSkyLines = !gsGetArgNamedBoolean(argc, argv, "--no-sky-lines"),
        .addSkyContinuum = !gsGetArgNamedBoolean(argc, argv, "--no-sky-continuum"),
        .addLunarContinuum = !gsGetArgNamedBoolean(argc, argv, "--no-lunar-continuum"),
        .addStrayLight = !gsGetArgNamedBoolean(argc, argv, "--no-stray-light"),
        .addDarkNoise = !gsGetArgNamedBoolean(argc, argv, "--no-dark"),
        .addReadoutNoise = !gsGetArgNamedBoolean(argc, argv, "--no-readout"),
    };
    return params;
}

int main(int argc, char* argv[]) {
    int i_arm;
    long i;
    double lambda, sample_factor;
    FILE* fp;

    SPECTRO_ATTRIB spectro;
    OBS_ATTRIB obs;
    
    double** sky;
    double** moon;
    double** stray;
    double** dark;
    double** readout;
    double** conv;
    double** atm_contop;
    double** atm_trans;
    double* temp;
  
    gsPrintCompilerFlags();

    /* Process command-line arguments */
    PARAMS params = getParams(argc, argv);

    /* Load config, do not degrade detector throughput */
    fp = gsOpenConfigFile(params.spectroConfig);
    gsReadSpectrographConfig(fp, &spectro, 1.0);
    gsCloseConfigFile(fp);

    fp = gsOpenConfigFile(params.obsConfig);
    gsReadObservationConfig(fp, &obs);
    gsCloseConfigFile(fp);
    
    // TODO: diffuse_stray should be in spectro config file but would need to change file format
    spectro.diffuse_stray = 0.02;

    /* Allocate noise vectors */
    gsAllocArmVectors(&spectro, &sky);
    gsAllocArmVectors(&spectro, &moon);
    gsAllocArmVectors(&spectro, &stray);
    gsAllocArmVectors(&spectro, &dark);
    gsAllocArmVectors(&spectro, &readout);
    gsAllocArmVectors(&spectro, &conv);
    gsAllocArmVectors(&spectro, &atm_contop);
    gsAllocArmVectors(&spectro, &atm_trans);

    /* Sky lines and continuum */
    if (params.addSkyLines) {
        printf("Computing sky lines ...\n");
        for (i_arm = 0; i_arm < spectro.N_arms; i_arm++) {
            printf(" //Arm%d//\n", i_arm);
            gsAddSkyLines(&spectro, &obs, i_arm, sky[i_arm]);
        }
    }

    if (params.addSkyContinuum) {
        printf("Computing sky continuum ...\n");
        for (i_arm = 0; i_arm < spectro.N_arms; i_arm++) {
            printf(" //Arm%d//\n", i_arm);
            gsAddSkyContinuum(&spectro, &obs, i_arm, sky[i_arm]);
        }
    }
    
    /* Moon */
    if (params.addLunarContinuum) {
        printf("Computing lunar continuum ...\n");
        for (i_arm = 0; i_arm < spectro.N_arms; i_arm++) {
            printf(" //Arm%d//\n", i_arm);
            gsAddLunarContinuum(&spectro, &obs, i_arm, moon[i_arm]);
        }
    }

    /* Stray light */
    /* Sum up sky and moon to compute stray light */
    if (params.addStrayLight) {
        printf("Computing stray light ...\n");
        temp = (double*)malloc((size_t)(MAXPIX*sizeof(double)));
        for (i_arm = 0; i_arm < spectro.N_arms; i_arm++) {
            sample_factor = gsGetSampleFactor(&spectro, i_arm);
            for(i = 0; i < spectro.npix[i_arm]; i++) {
                temp[i] = (sky[i_arm][i] + moon[i_arm][i]) / sample_factor;
            }
            gsAddStrayLight(&spectro, &obs, i_arm, stray[i_arm], temp);
        }
        free(temp);
    }

    /* Dark and read-out noise */
    if (params.addDarkNoise) {
        printf("Computing dark and readout noise ...\n");
        for (i_arm = 0; i_arm < spectro.N_arms; i_arm++) {
            printf(" //Arm%d//\n", i_arm);
            gsAddDarkNoise(&spectro, &obs, i_arm, dark[i_arm]);
            gsAddReadoutNoise(&spectro, i_arm, readout[i_arm]);
        }
    }

    /* Compute overall transmission function including reddening, atmosphere and instrument */
    /* Do conversion between erg/cm2/s to counts per pixel */
    printf("Computing transmission functions ...\n");
    for (i_arm = 0; i_arm < spectro.N_arms; i_arm++) {
        printf(" //Arm%d//\n", i_arm);
        for(i = 0; i < spectro.npix[i_arm]; i++) {
            lambda = spectro.lmin[i_arm] + spectro.dl[i_arm] * i;
            atm_contop[i_arm][i] = gsAtmContOp(&obs, lambda);
            atm_trans[i_arm][i] = gsAtmTransInst(&spectro, &obs, i_arm, lambda);
            conv[i_arm][i] = gsConversionFunction(&spectro, &obs, i_arm, lambda);
        }
    }
    
    /* Write into output file */
    fp = gsOpenOutputFile(params.noiseFile);
    gsWriteObservationConfig(fp, &obs, "# ");
    fprintf(fp, "# arm i wave sky moon stray dark readout conv sample_factor atm_contop atm_trans \n");
    for (i_arm = 0; i_arm < spectro.N_arms; i_arm++) {
        sample_factor = gsGetSampleFactor(&spectro, i_arm);
        for (i = 0; i < spectro.npix[i_arm]; i++) {
            lambda = spectro.lmin[i_arm] + spectro.dl[i_arm] * (i + 0.5);
	        fprintf(fp, "%1d %4ld %7.4lf %11.5le %11.5le %11.5le %11.5le %11.5le %11.5le %11.5le %11.5le %11.5le\n", 
                spectro_arm(&spectro, i_arm), i, lambda, 
                sky[i_arm][i], moon[i_arm][i], stray[i_arm][i], dark[i_arm][i], readout[i_arm][i], 
                conv[i_arm][i], sample_factor, atm_contop[i_arm][i], atm_trans[i_arm][i]);
      }
      fprintf(fp, "\n");
    }
    gsCloseOutputFile(fp);
    printf(" --> Done.\n");

    /* De-allocate noise vectors */
    gsFreeArmVectors(&spectro, sky);
    gsAllocArmVectors(&spectro, &moon);
    gsAllocArmVectors(&spectro, &stray);
    gsAllocArmVectors(&spectro, &dark);
    gsAllocArmVectors(&spectro, &readout);
    gsAllocArmVectors(&spectro, &conv);
    gsAllocArmVectors(&spectro, &atm_contop);
    gsAllocArmVectors(&spectro, &atm_trans);
}