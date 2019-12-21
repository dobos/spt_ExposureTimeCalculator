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
} PARAMS;

PARAMS getParams(int argc, char* argv[]) {
    PARAMS params;
    params.spectroConfig = argv[1];
    params.obsConfig = argv[2];
    params.noiseFile = argv[3];
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

    /* Sky lines and continuum */
    printf("Computing sky lines and continuum ...\n");
    for (i_arm = 0; i_arm < spectro.N_arms; i_arm++) {
        printf(" //Arm%d//\n", i_arm);
        gsAddSkyLines(&spectro, &obs, i_arm, sky[i_arm]);
        gsAddSkyContinuum(&spectro, &obs, i_arm, sky[i_arm]);
    }

    /* Moon */
    printf("Computing lunar continuum ...\n");
    for (i_arm = 0; i_arm < spectro.N_arms; i_arm++) {
        printf(" //Arm%d//\n", i_arm);
        gsAddLunarContinuum(&spectro, &obs, i_arm, moon[i_arm]);
    }

    /* Stray light */
    /* Sum up sky and moon to compute stray light */
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

    /* Dark and read-out noise */
    printf("Computing dark and readout noise ...\n");
    for (i_arm = 0; i_arm < spectro.N_arms; i_arm++) {
        printf(" //Arm%d//\n", i_arm);
        gsAddDarkNoise(&spectro, &obs, i_arm, dark[i_arm]);
        gsAddReadoutNoise(&spectro, i_arm, readout[i_arm]);
    }

    /* Compute overall transmission function including reddening, atmosphere and instrument */
    /* Do conversion between erg/cm2/s to counts per pixel */
    printf("Computing transmission function ...\n");
    for (i_arm = 0; i_arm < spectro.N_arms; i_arm++) {
        printf(" //Arm%d//\n", i_arm);
        for(i = 0; i < spectro.npix[i_arm]; i++) {
            lambda = spectro.lmin[i_arm] + spectro.dl[i_arm] * i;
            conv[i_arm][i] = gsConversionFunction(&spectro, &obs, i_arm, lambda);
        }
    }
    
    /* Write into output file */
    fp = gsOpenOutputFile(params.noiseFile);
    gsWriteObservationConfig(fp, &obs, "# ");
    fprintf(fp, "# arm i wave sky moon stray dark readout conv sample_factor\n");
    for (i_arm = 0; i_arm < spectro.N_arms; i_arm++) {
        sample_factor = gsGetSampleFactor(&spectro, i_arm);
        for (i = 0; i < spectro.npix[i_arm]; i++) {
            lambda = spectro.lmin[i_arm] + spectro.dl[i_arm] * (i + 0.5);
	        fprintf(fp, "%1d %4ld %7.4lf %11.5le %11.5le %11.5le %11.5le %11.5le %11.5le %11.5le\n", 
                spectro_arm(&spectro, i_arm), i, lambda, 
                sky[i_arm][i], moon[i_arm][i], stray[i_arm][i], dark[i_arm][i], readout[i_arm][i], 
                conv[i_arm][i], sample_factor);
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
}