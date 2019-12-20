#ifndef __LIBGSETC_H__
#define __LIBGSETC_H__

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/* VERSION 5 */
/* MODIFICATION 6 */
/* [OII] recovery histogram data */
#define ZMIN_OII 0.10
#define NZ_OII 24
#define DZ_OII 0.10

/* Conversions */
#define ARCSEC_PER_URAD 0.206264806247097
#define DEGREE 0.017453292519943278
#define PHOTONS_PER_ERG_1NM 5.03411747e8

/* Ratio of half-light to scale radius for exponential profile galaxy */
#define RAT_HL_SL_EXP 1.67834

/* Maximum legal number of arms */
#define MAXARM 3

/* Maximum number of l-pixels per spectrograph arm */
#define MAXPIX 8192

/* Maximum length of throughput table */
#define MAXNTHR 1024

/* Spectrograph PSF length (must be even) */
#define SP_PSF_LEN 64

/* Spectrograph attributes structure */
typedef struct {
  double D_outer, EFL[5];     /* Outer diameter & effective focal length (meters) */
  double fiber_ent_rad;       /* Fiber entrance radius (microns) */
  double centobs;             /* Central obscuration */
  double rfov;                /* Field of view radius (degrees) */
  double rms_spot[5];         /* RMS spot size per axis: 5 values from center to edge (microns) */
  double vignette[5];         /* Vignetting factor: 5 values from center to edge */
  int N_arms;                 /* Number of arms */
  int MR;                     /* True iff we are configured to use the medium resolution grating */

  double lmin[MAXARM];        /* Min wavelength in nm */
  double lmax[MAXARM];        /* Max wavelength in nm */
  long npix[MAXARM];          /* Number of pix from min to max wavelength */
  double dl[MAXARM];          /* Wavelength spacing, nm/pix, derived */
  long width[MAXARM];         /* Width of trace used in actual analysis, pixels */
  double fratio[MAXARM];      /* Camera f/ratio in ith arm */
  double thick[MAXARM];       /* Camera detector thickness (microns) */
  double pix[MAXARM];         /* Camera pixel scale (microns) */
  double temperature[MAXARM]; /* Camera temperature (K) */
  double rms_cam[MAXARM];     /* Camera rms spot per axis (microns) */
  double diam[MAXARM];        /* Geometric fiber image diameter (microns) */
  double dark[MAXARM];        /* Dark current (e/pix/s) */
  double read[MAXARM];        /* Read noise (e rms) */
  double sep[MAXARM];         /* Trace spacing (microns) */
  double nline[MAXARM];       /* Effective number of lines */

  int N_thr;                  /* Number of points in throughput grid */
  int istart[MAXARM+1];       /* Position at which the i_arm grid begins; istart[N_arms] = N_thr */
  double l[MAXNTHR];          /* Wavelength at grid points (nm) */
  double T[MAXNTHR];          /* Throughput at grid points (excluding atmosphere & fiber geom) */

  int Dtype[MAXARM];          /* =0 for Si, 1 for HgCdTe */

  double sysfrac;             /* Sky subtraction systematic (per pixel per exposure) */
  double diffuse_stray;       /* Diffuse stray light amplitude */

} SPECTRO_ATTRIB;

static int spectro_arm(const SPECTRO_ATTRIB *spectro, int ia)
{
   if (!spectro->MR) {
      return ia;
   } else {
      return (ia == 1) ? 3 : ia;
   }
}

/* Observing condition attributes structure */
typedef struct {
  double seeing_fwhm_800; /* Seeing FWHM @ 800 nm */
  double elevation;       /* Elevation in meters above sea level */
  double zenithangle;     /* Angle from the zenith, in degrees */
  double lunarphase;      /* Lunar phase: 0.0 (new Moon) --> 0.5 (full Moon) --> 1.0 (new Moon) */
  double lunarangle;      /* Angle from the Moon to the line of sight, in degrees */
  double lunarZA;         /* Angle from the Moon to the zenith, in degrees; set to >90 for no Moon */
  double EBV;             /* Dust reddening column, E(B-V) */
  unsigned long skytype;  /* Bitmask for assumptions on atmospheric conditions */
} OBS_ATTRIB;

/* Added by Y.Moritani for input mag. file: 20150422 */
/* Modified by L.Dobos to use struct and pass around variable> 20191219 */
typedef struct {
  double* lambda_inmag;
  double* mag_inmag;
  int num_inmag;
} MAGFILE;
/* Added by Y.Moritani for input mag. file: 20150422 : end*/

/* Exported functions */

double gsGeometricThroughput(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, double lambda,
  double r_eff, double decent, double fieldang, unsigned long flags);

double gsAeff(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double lambda, double fieldang);

double gsFracTrace(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double lambda, int tr);

double gsAtmTransInst(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double lambda, unsigned long flags);

double gsConversionFunction(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, 
  double r_eff, double decent, double fieldang, double t_exp, double lambda, unsigned long flags);

double gsMagToFlux(double mag);
double gsEBVToAttn(OBS_ATTRIB* obs, double lambda);

double gsGetSampleFactor(SPECTRO_ATTRIB* spectro, int i_arm);

void gsGetNoise(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double fieldang,
  double *Noise, double *SkyMod, double t_exp, unsigned long flags);

double gsGetSNR_OII(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double z,
  double F, double sigma_v, double r_eff, double src_cont, double ROII, double decent,
  double fieldang, double *Noise, double t_exp, unsigned long flags, int snrType);

double gsGetSNR_Single(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double mag, double lambda,
  double F, double sigma_v, double r_eff, double decent, double fieldang, double *Noise,
  double t_exp, unsigned long flags, int snrType, MAGFILE* magfile2);

void gsGetSNR_Continuum(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double mag,
  double r_eff, double decent, double fieldang, double *Noise, double t_exp, unsigned long flags,
  MAGFILE* magfile2,
  double *out_SNR_curve, double *out_count_curve, double *out_noise_curve, double *out_mag_curve, double *out_trans_curve, double *out_sample_factor_curve);

void gsAddSkyLines(SPECTRO_ATTRIB* spectro, OBS_ATTRIB* obs, int i_arm, double* Noise,
  double fieldang, double t_exp, unsigned long flags);

void gsAddLunarContinuum(SPECTRO_ATTRIB* spectro, OBS_ATTRIB* obs, int i_arm,
  double* Noise, double fieldang, double t_exp, unsigned long flags);

void gsAddSkyContinuum(SPECTRO_ATTRIB* spectro, OBS_ATTRIB* obs, int i_arm,
  double* Noise, double fieldang, double t_exp, unsigned long flags);

void gsAddStrayLight(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm,
  double* Noise, double* sky);

void gsAddDarkNoise(SPECTRO_ATTRIB *spectro, int i_arm, double t_exp, double* Noise);

void gsAddReadoutNoise(SPECTRO_ATTRIB *spectro, int i_arm, double* Noise);

void gsAllocArmVectors(SPECTRO_ATTRIB* spectro, double*** spec);
void gsFreeArmVectors(SPECTRO_ATTRIB* spectro, double** spec);
void gsReadMagfile(MAGFILE* magfile, char* filename);
void gsCopyMagfile(MAGFILE* magfile, MAGFILE* magfile2);
double gsInterpolateMagfile(MAGFILE* magfile, double lambda);
void gsPrintCompilerFlags();
FILE* gsOpenConfigFile(const char* filename);
void gsCloseConfigFile(FILE* fp);
FILE* gsOpenOutputFile(const char* filename);
void gsCloseOutputFile(FILE* fp);
void gsReadSpectrographConfig(FILE*, SPECTRO_ATTRIB *spectro, double degrade);
void gsReadObservationConfig(FILE *fp, OBS_ATTRIB *obs, double* fieldang, double* t_exp);
void gsReadObsConfig_Legacy(OBS_ATTRIB* obs, SPECTRO_ATTRIB* spectro, double* fieldang, double* decent, double* t, int* n_exp);

#endif