#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "libgsetc.h"

/* Generates the signal/noise ratio for a single emission line given the noise vector Noise taking into consideration the continuum effect.
 *
 * Line parameters are lambda (nm), F (erg/cm2/s), sigma_v (km/s), r_eff (arcsec),
 * decent (arcsec), fieldang (degrees).
 *
 * The "snrTypes" types are:
 *  0 = 1D optimal
 *  1 = uniform matched filter
 */
double gsGetSNR_Single(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double mag, double lambda,
  double F, double sigma_v, double *Noise, int snrType, MAGFILE* magfile2) {

  long ipix,Npix;
  double SNR = 0;
  double numer, denom;
  double *Signal;
  /* Added by K.Yabe for input mag. file: 20160205 */
  double counts, trans, den, x;
  int flag = 0;
  double src_cont;

  if(mag==-99.9) flag=1;
  /* Added by K.Yabe for input mag. file: 20160205 : end */

  /* Added by Y.Moritani for input mag. file: 20150422 */
  /* interpolate magnitude for a given lambda */
  /* mag is used as a fllag: -99.9 ... use input file*/
  if (flag) {
    mag = gsInterpolateMagfile(magfile2, lambda);
  }

  /* Atmospheric transmission */
  trans = den = 0.;
  for(x=-4;x<4.01;x+=.2) {
    trans += gsAtmTrans(obs,lambda) * exp(-0.5*x*x);
    den += exp(-0.5*x*x);
  }
  trans /= den;
 
  /* Determine how many counts per Hz we get from the object continuum */
  src_cont = 3.631e-20 * pow(10., -0.4*mag); /* in erg/cm2/s */
  counts = src_cont * trans
             * pow(10., -0.4*gsGalactic_Alambda__EBV(lambda)*obs->EBV)
             * gsGeometricThroughput(spectro,obs,lambda)
             * gsFracTrace(spectro,i_arm,lambda,0)
             * PHOTONS_PER_ERG_1NM * lambda * obs->t_exp 
             * gsAeff(spectro,obs,i_arm) 
             * gsThroughput(spectro,i_arm,lambda) * 1e4;
  #ifdef HGCDTE_SUTR
    if (spectro->Dtype[i_arm]==1) counts *= 1.2;
  #endif
  /* Convert from per Hz --> l-per pixel */
  counts *= 2.99792458e17*spectro->dl[i_arm]/(lambda*lambda);

  /* Allocate and get the signal vector */
  Npix = spectro->npix[i_arm];
  Signal = (double*)malloc((size_t)(Npix*sizeof(double)));
  gsGetSignal(spectro,obs,i_arm,lambda,F,sigma_v,Signal);

  /* 1D optimal */
  if (snrType == 0) {
    SNR = 0;
    for(ipix=0;ipix<Npix;ipix++) SNR += Signal[ipix]*Signal[ipix]/(counts+Noise[ipix]);
    SNR = sqrt(SNR);
  }

  /* uniform matched filter */
  if (snrType == 1) {
    numer = denom = 0;
    for(ipix=0;ipix<Npix;ipix++) {
      numer += Signal[ipix]*Signal[ipix];
      denom += Signal[ipix]*Signal[ipix]*(counts+Noise[ipix]);
    }
    SNR = numer>=0.001? numer/sqrt(denom): 0;
  }

  free((char*)Signal);
  return(SNR);
}

/* Obtains the SNR for the [OII] doublet.
 * Line parameters are z, F (erg/cm2/s), sigma_v (km/s), r_eff (arcsec),
 * src_cont (erg/cm2/s/Hz), and ROII.
 * Observing parameters are decent (arcsec), fieldang (degrees).
 *
 * The "snrTypes" types are:
 *  0 = 1D optimal, brighter feature, assuming 1:1 ratio
 *  1 = uniform matched filter, brighter feature, assuming 1:1 ratio
 */
/*
 * in main loop
 * snr[ia] = gsGetSNR_OII(&spectro,&obs,ia,z,1e-16,70.,REF_SIZE,0.,1.,decent,fieldang,spNoise[ia],t,0x0,snrType)
 * *sqrt((double)n_exp);
 */
double gsGetSNR_OII(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double z,
  double F, double sigma_v, double src_cont, double ROII, double *Noise, int snrType) {

  int i;
  double SNR=0.;
  double lambda[2], indivSNR[2], frac[2];
  long Npix, ipix;
  double *Signal0, *Signal1, *myNoise;
#ifndef NO_OBJ_CONTINUUM
  double counts, ll, trans, den, x;
#endif

  lambda[0] = 372.71*(1.+z);
  lambda[1] = 372.98*(1.+z);

  /* in this simulation ROII = 1. ; comment by Y. Moritani*/
  if (ROII<0.667) ROII=0.667;
  if (ROII>3.87) ROII=3.87;

  frac[0] = ROII/(1.+ROII);
  frac[1] = 1./(1.+ROII);

  /* Build noise vector including continuum */
  Npix = spectro->npix[i_arm];
  myNoise = (double*)malloc((size_t)(Npix*sizeof(double)));
  for(ipix=0;ipix<Npix;ipix++) myNoise[ipix]=Noise[ipix];
#ifndef NO_OBJ_CONTINUUM
  /* Atmospheric transmission */
  trans = den = 0.;
  for(x=-4;x<4.01;x+=.2) {
    trans += gsAtmTrans(obs,lambda[0] + (0.5+0.5*x)*(lambda[1]-lambda[0])) * exp(-0.5*x*x);
    den += exp(-0.5*x*x);
  }
  trans /= den;
 
  /* Determine how many counts per Hz we get from the object continuum */
  ll = (lambda[0]+lambda[1])/2.;
  counts = src_cont * trans
           * pow(10., -0.4*gsGalactic_Alambda__EBV(ll)*obs->EBV)
           * gsGeometricThroughput(spectro,obs,ll)
           * gsFracTrace(spectro,i_arm,ll,0)
           * PHOTONS_PER_ERG_1NM * ll * obs->t_exp 
           * gsAeff(spectro,obs,i_arm)
           * gsThroughput(spectro,i_arm,ll) * 1e4;
#ifdef HGCDTE_SUTR
  if (spectro->Dtype[i_arm]==1) counts *= 1.2;
#endif
  /* Convert from per Hz --> l-per pixel */
  counts *= 2.99792458e17*spectro->dl[i_arm]/(ll*ll);

  /* Add to all points in the data vector ... not really correct everywhere but good
   * near the line where it matters.
   */
  for(ipix=0;ipix<Npix;ipix++) myNoise[ipix] += counts;
#endif

  /* 1D optimal - brighter feature */
  if (snrType == 0) {
    for(i=0;i<2;i++)
      indivSNR[i] = gsGetSNR(spectro,obs,i_arm,lambda[i],frac[i]*F,sigma_v,myNoise,0);
    SNR = indivSNR[0]>indivSNR[1]? indivSNR[0]: indivSNR[1];
  }

  /* uniform matched filter - brighter feature */
  if (snrType == 1) {
    for(i=0;i<2;i++)
      indivSNR[i] = gsGetSNR(spectro,obs,i_arm,lambda[i],frac[i]*F,sigma_v,myNoise,1);
    SNR = indivSNR[0]>indivSNR[1]? indivSNR[0]: indivSNR[1];
  }

  /* Combined "optimal" of the 2 lines */
  if (snrType == 2) {
    Npix = spectro->npix[i_arm];
    Signal0 = (double*)malloc((size_t)(Npix*sizeof(double)));
    Signal1 = (double*)malloc((size_t)(Npix*sizeof(double)));
    gsGetSignal(spectro,obs,i_arm,lambda[0],frac[0]*F,sigma_v,Signal0);
    gsGetSignal(spectro,obs,i_arm,lambda[1],frac[1]*F,sigma_v,Signal1);
    SNR = 0;
    for(ipix=0;ipix<Npix;ipix++) SNR += (Signal0[ipix]+Signal1[ipix])*(Signal0[ipix]+Signal1[ipix])/myNoise[ipix];
    SNR = sqrt(SNR);
    free((char*)Signal0);
    free((char*)Signal1);
  }

  free((char*)myNoise);
  return(SNR);
}

/* Obtains the continuum S/N per pixel per exposure as a function of the source magnitude (AB) in the
 * specified arm. Output is to out_SNR_curve[spectro->npix[i_arm]].
 */

/* Modified by Y.Moritani for input mag. file: 20150422 :*/
void gsGetSNR_Continuum(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double mag,
  double *Noise, MAGFILE* magfile2,
  double *out_SNR_curve, double *out_count_curve, double *out_noise_curve, double *out_mag_curve, double *out_trans_curve, double *out_sample_factor_curve) {

  long ipix;
  double lambda, src_cont, counts;
  double sample_factor;

  /* Added by Y.Moritani for input mag. file: 20150422 */
  int flag = 0;

  if(mag==-99.9) flag=1;
  /* Added by Y.Moritani for input mag. file: 20150422 : end */

  /* Pre-factor for sampling the Poisson distribution */
  sample_factor = 1.0;
#ifdef HGCDTE_SUTR
  if (spectro->Dtype[i_arm]==1) sample_factor = 1.2;
#endif

  for(ipix=0;ipix<spectro->npix[i_arm];ipix++) {
    gsPrintProgress(spectro->npix[i_arm], ipix);
    lambda = spectro->lmin[i_arm] + spectro->dl[i_arm]*ipix;
 
    /* Added by Y.Moritani for input mag. file: 20150422 */
    /* interpolate magnitude for a given lambda */
    /* mag is used as a fllag: -99.9 ... use input file*/
    if(flag) {  
      mag = gsInterpolateMagfile(magfile2, lambda);
    }
    
    /* Determine how many counts per Hz we get from the object continuum */
    src_cont = gsMagToFlux(mag); /* in erg/cm2/s */
    counts = src_cont * gsConversionFunction(spectro,obs,i_arm,lambda);

    /* Report S/N ratio */
    out_SNR_curve[ipix] = counts/sqrt(sample_factor*counts + Noise[ipix]);
    /* Modified by K. Yabe 20150525 */
    out_count_curve[ipix] = counts;
    out_noise_curve[ipix] = sample_factor*counts + Noise[ipix];
    out_mag_curve[ipix] = mag;
    out_trans_curve[ipix] = counts / src_cont;
    out_sample_factor_curve[ipix] = sample_factor;
  }
  return;
}

/* === MAIN PROGRAM === */

/* Reads observation parameters from stdin.
 * Called by the python wrapper
 */
int main(void) {
  FILE *fp, *fq;
  int ia;
  long i, j, k;
  SPECTRO_ATTRIB spectro;
  OBS_ATTRIB obs;
  char FileName[256], OutFileNoise[256], OutFileSNR[256], OutFileSNRAB[256];
  char InFileOII[256], OutFileOII[256];
  int flag_reused;
  double lambda, z;
  double snr[MAXARM], snrtot, Aeff;
  double **spNoise;
  double **spSky;
  long id;
  double ROII, FOII, contOII, sigma;
  double min_SNR = 0.;
  int snrType;
  long ngal[NZ_OII], ngtot;
  double snrmax16, mdlf;
  double *snrcont;
  double *snrcontcount;
  double *snrcontnoise;
  double *magcont;
  double *snctrans;
  double *samplefac;
  int arm;
  long pix;
  double wav, innoise, inskymod;
  int proc, proc_tot;
  char buf[256];

  MAGFILE inmag, inmag2;

  /* Added by Y.Moritani for input mag. file: 20150422*/
  char InFileMag[256];
  double mag = 22.5;  /* default value of mag, used as a flag for inputfile*/

  /* Added by Y.Moritani for line SNR. file: 20150427*/
  char OutFileSNR2[256];
  double flux_emi, sigma_emi;
  double degrade;

  /* Tell us what flags are on */
  gsPrintCompilerFlags();

  /* Intializations */
  for(j=0;j<NZ_OII;j++) ngal[j] = 0;

  /* Get the spectrograph properties */
  //printf("Enter spectrograph configuration file: ");
  if (!scanf("%255s", FileName)) {
    fprintf(stderr, "Error: Can't get input file.\n");
    exit(1);
  }

  if(scanf("%lg", &(degrade))==EOF)
    degrade = 1.0;

  fp = gsOpenConfigFile(FileName);
  gsReadSpectrographConfig(fp, &spectro, degrade);
  gsCloseConfigFile(fp);

  /* Allocate noise vectors */
  gsAllocArmVectors(&spectro, &spNoise);
  gsAllocArmVectors(&spectro, &spSky);

  gsReadObsConfig_Legacy(&obs, &spectro);

  /* Output files */
  //printf("Noise data reused?: [1=yes/0=no] ");
  if(scanf("%d", &flag_reused)==EOF)
    flag_reused=0;
  //printf("Enter output file for noise vector: ");
  if (!scanf("%255s", OutFileNoise)) {
    fprintf(stderr, "Error: Failed to read name of output file.\n");
    exit(1); 
  }  
  //printf("Enter output file for ELG S/N curve: [- for no ELG S/N curve output] ");
  if (!scanf("%255s", OutFileSNR)) {
    fprintf(stderr, "Error: Failed to read name of output file.\n");
    exit(1); 
  }

  /* Added by Y.Moritani for line SNR. file: 20150427 */
  //printf("Enter output file for ELG S/N curve2: [- for no ELG S/N curve2 output] ");
  if (!scanf("%255s", OutFileSNR2)) {
    fprintf(stderr, "Error: Failed to read name of output file.\n");
    exit(1); 
  } 
  //printf("Enter flux of the emission line [erg cm-2 s-1]: ");
  if(scanf("%lg", &flux_emi)==EOF)
    flux_emi=1.0e-7;
  //printf("Enter velocity width of the emission line [km s-1]: ");
  if(scanf("%lg", &sigma_emi)==EOF)
    sigma_emi=70.;
  /* Added by Y.Moritani for line SNR. file: 20150427 : end */

  //printf("Enter output file for continuum curve: [- for no continuum S/N curve output] ");
  if (!scanf("%255s", OutFileSNRAB)) {
    fprintf(stderr, "Error: Failed to read name of output file.\n");
    exit(1); 
  }

  /* Files for [OII] detection */
  //printf("Enter [OII] input catalogue file: [- for no OII computation] ");
  if (!scanf("%255s", InFileOII)) {
    fprintf(stderr, "Error: Failed to read name of input file.\n");
    exit(1); 
  }  
  if (!scanf("%255s", OutFileOII)) {
      fprintf(stderr, "Error: Failed to read name of output file.\n");
      exit(1); 
  }
    //printf("Enter minimum SNR: ");
  if (scanf("%lg", &min_SNR)==EOF)  {
      min_SNR=9.;
  }

  /* Added by Y.Moritani for input mag. file: 20150422*/
  /* File for Input mag:*/
  //printf("Enter magnitude input file: [- for no designated file] ");
  if (!scanf("%255s", InFileMag)) {
    fprintf(stderr, "Error: Failed to read name of input file.\n");
    exit(1); 
  }
  if (strcmp("-",InFileMag)!=0) {
    gsReadMagfile(&inmag, InFileMag);
    gsCopyMagfile(&inmag, &inmag2);
    mag=-99.9; /* used as a flag for inputfile*/
  }
  /* Added by Y.Moritani for input mag. file: 20150422 : end*/
  //printf("Enter the effective radius of the galaxy [arcsec]:\n");
  if(scanf("%lg", &obs.r_eff)==EOF)
    obs.r_eff=0.3;
  /* Added by K.Yabe */

  /* Encircled energy in fiber */
  printf("Fiber aperture factor [@800nm, r_eff=%.2lf\"(exp)] = %10.8lf\n", obs.r_eff, gsGeometricThroughput(&spectro, &obs, 800));
  printf("Fiber aperture factor [@800nm,     point source] = %10.8lf\n", gsGeometricThroughput(&spectro, &obs, 800));

  proc=0;
  proc_tot=0;
  if (flag_reused == 0) {
    proc_tot+=1;
  }
  if (strcmp("-",OutFileSNR)!=0) {
    proc_tot+=1;
  }
  if (strcmp("-",OutFileSNR2)!=0) {
    proc_tot+=1;
  }
  if (strcmp("-",OutFileSNRAB)!=0) {
    proc_tot+=1;
  }
  if (strcmp("-",InFileOII)!=0) {
    proc_tot+=1;
  }
  /* Generate and write noise vector */
  if (flag_reused == 1) {
    printf("Loading noise vector ...\n");
    k=0;
    fp=fopen(OutFileNoise,"r");
    while(fgets(buf, 256, fp) != NULL){
      if(strncmp(buf,"\n",1) == 0) continue;
      sscanf(buf, "%d %ld %lf %le %le", &arm, &pix, &wav, &innoise, &inskymod);
      if (k>=spectro.npix[0]+spectro.npix[1]) {
        spNoise[2][k-spectro.npix[0]-spectro.npix[1]]=innoise;
        spSky[2][k-spectro.npix[0]-spectro.npix[1]]=inskymod;
      }
      else if (k>=spectro.npix[1]) {
        spNoise[1][k-spectro.npix[0]]=innoise;
        spSky[1][k-spectro.npix[0]]=inskymod;
      }
      else {
        spNoise[0][k]=innoise;
        spSky[0][k]=inskymod;
      }
      k++;
    }
  }
  else {
    proc+=1;
    printf("(%d/%d) Computing noise vector ...\n",proc, proc_tot);
    for(ia=0;ia<spectro.N_arms;ia++)
      gsGetNoise(&spectro,&obs,ia,spNoise[ia],spSky[ia]);
    fp = fopen(OutFileNoise, "w");
    for(ia=0;ia<spectro.N_arms;ia++) {
      for(i=0;i<spectro.npix[ia];i++) {
	      fprintf(fp, "%1d %4ld %7.4lf %11.5le %11.5le\n", spectro_arm(&spectro, ia),
		      i, lambda=spectro.lmin[ia]+spectro.dl[ia]*(i+0.5), spNoise[ia][i], spSky[ia][i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
    printf(" --> Done.\n");
  }

  /* Generate and write the S/N ratio vector for the [OII] doublet as a function of redshift,
   * at a fiducial case of r_eff=0.3", F=1e-16 erg/cm2/s, and sigma_v=70 km/s.
   */

  snrType = 0;

#define REF_SIZE 0.30
  if (strcmp("-",OutFileSNR)!=0) {
    proc+=1;
    snrType = 2;
    printf("(%d/%d) Computing SNR curve for [OII] emission lines ...\n",proc,proc_tot);
    fp = fopen(OutFileSNR, "w");
    for(z=0.1;z<2.5001;z+=0.0002) {
      printf("      --> %.0f percent done ...\r",41.666*(z-0.1));
      snrtot = 0.;
      for(ia=0;ia<spectro.N_arms;ia++) {
        snr[ia] = 0.;
        if (spectro.lmin[ia]<373.8*(1+z) && 371.8*(1+z)<spectro.lmax[ia])
          snr[ia] = gsGetSNR_OII(&spectro,&obs,ia,z,flux_emi,sigma_emi,0.,1.,spNoise[ia],snrType)
                    *sqrt((double)obs.n_exp);
        snrtot += snr[ia]*snr[ia];
      }
      snrtot = sqrt(snrtot);

      Aeff = 0.;
      for(ia=0;ia<spectro.N_arms;ia++)
        if (spectro.lmin[ia]<372.845*(1+z) && 372.845*(1+z)<spectro.lmax[ia])
          Aeff += gsAeff(&spectro,&obs,ia)
                  * gsThroughput(&spectro,ia,372.845*(1+z));

      fprintf(fp, "%6.4lf %7.2lf %7.2lf %8.6lf %8.5lf",
        z, 372.71*(1+z), 372.98*(1+z), gsGeometricThroughput(&spectro, &obs, 372.845*(1+z)), Aeff
      );
      for(ia=0;ia<spectro.N_arms;ia++)
        fprintf(fp, " %8.4lf", snr[ia]);
      fprintf(fp, " %8.4lf\n", snrtot);
    }
    printf("\n");
    printf(" --> Done.\n");
    fclose(fp);
  }

  /* Added by Y.Moritani for line SNR. file: 20150427 */
  /* Modified by K.Yabe for line SNR. file: 20160126 */
  if (strcmp("-",OutFileSNR2)!=0) {
    proc+=1;
    snrType = 0;
    printf("(%d/%d) Computing SNR curve for a single line with f=%.2e [erg cm-2 s-1], sigma=%.0lf [km s-1] ...\n", proc, proc_tot, flux_emi, sigma_emi);
    fp = fopen(OutFileSNR2, "w");
    for(z=0.1;z<2.7627;z+=0.0002) {
      printf("      --> %.0f percent done ...\r",37.556*(z-0.1));
      snrtot = 0.;
      for(ia=0;ia<spectro.N_arms;ia++) {
        snr[ia] = 0.;
        if (spectro.lmin[ia]<345.5*(1+z) && 345.5*(1+z)<spectro.lmax[ia])
          snr[ia] = gsGetSNR_Single(&spectro,&obs,ia,mag,345.5*(1+z),
          flux_emi,sigma_emi,spNoise[ia],0,&inmag2)*sqrt((double)obs.n_exp);
        snrtot += snr[ia]*snr[ia];
      }
      snrtot = sqrt(snrtot);

      Aeff = 0.;
      for(ia=0;ia<spectro.N_arms;ia++)
        if (spectro.lmin[ia]<345.5*(1+z) && 345.5*(1+z)<spectro.lmax[ia])
          Aeff += gsAeff(&spectro,&obs,ia)
                  * gsThroughput(&spectro,ia,345.5*(1+z));

/*      fprintf(fp, "%6.4lf %7.2lf %8.6lf %8.5lf", z, 345.5*(1+z), gsGeometricThroughput(&spectro, &obs, 345.5*(1+z), ref_input, decent, 0, 0x0), Aeff); */
      fprintf(fp, "%7.2lf %8.6lf %8.5lf", 345.5*(1+z), gsGeometricThroughput(&spectro, &obs, 345.5*(1+z)), Aeff);
      for(ia=0;ia<spectro.N_arms;ia++)
        fprintf(fp, " %8.4lf", snr[ia]);
      fprintf(fp, " %8.4lf\n", snrtot);
    }
    printf("\n");
    printf(" --> Done.\n");
    fclose(fp);
  }
  /* Added by Y.Moritani for line SNR. file: 20150427 : end */
#undef REF_SIZE

  /* Generate and write the continuum S/N for a magnitude 22.5 AB object.
   */
  if (strcmp("-",OutFileSNRAB)!=0) {
    proc+=1;
    snrType = 0;
    printf("(%d/%d) Computing SNR curve for continuum ...\n",proc,proc_tot);
    fp = fopen(OutFileSNRAB, "w");
    for(ia=0;ia<spectro.N_arms;ia++) {
      /* Allocate data vectors for the current arm */
      snrcont = (double*)malloc((size_t)(spectro.npix[ia]));
      snrcontcount = (double*)malloc((size_t)(spectro.npix[ia]));
      snrcontnoise = (double*)malloc((size_t)(spectro.npix[ia]));
      magcont = (double*)malloc((size_t)(spectro.npix[ia]));
      snctrans = (double*)malloc((size_t)(spectro.npix[ia]));
      samplefac = (double*)malloc((size_t)(spectro.npix[ia]));

      /* Modified by Y. Moritani for input mag. file: 20150422 :*/
      /* Modified by K. Yabe for counts output: 20150525 :*/
      //gsGetSNR_Continuum(&spectro,&obs,ia,22.5,0.0,decent,fieldang,spNoise[ia],t,0x0,snrcont);
      gsGetSNR_Continuum(&spectro,&obs,ia,mag,spNoise[ia],&inmag2,
        snrcont,snrcontcount,snrcontnoise,magcont,snctrans,samplefac);
      for(j=0;j<spectro.npix[ia];j++) {
        fprintf(fp, "%1d %4ld %9.3lf %8.4lf %11.5lE %11.5lE %11.5lE %11.5lE %11.5lE %11.5lE  %11.5lE\n",
		            spectro_arm(&spectro, ia), j, spectro.lmin[ia]+spectro.dl[ia]*j,snrcont[j]*sqrt((double)obs.n_exp),snrcontcount[j],spNoise[ia][j],snrcontnoise[j],magcont[j],snctrans[j],samplefac[j],spSky[ia][j]);
      }

      free(snrcont);
      free(snrcontcount);
      free(snrcontnoise);
      free(magcont);
      free(snctrans);
      free(samplefac);
    }
    printf("\n");
    printf(" --> Done.\n");
    fclose(fp);
}

  /* [OII] catalog file */
  if (strcmp("-",InFileOII)!=0) {
    proc+=1;
    snrType = 2;
    printf("(%d/%d) Processing [OII] emitter catalog ...\n", proc, proc_tot);

    /* Determine MDLF for point sources -- use this later to eliminate from the catalog
     * all objects where there is simply no way they could be detected.
     * Limits are computed for 1:1 line ratio, and reduced by 1.6 since neither line can have
     * more than 80% of the doublet flux due to atomic physics considerations.
     */
    snrmax16 = 1;
    for(z=0.1;z<2.5001;z+=0.00025) {
      snrtot = 0.;
      for(ia=0;ia<spectro.N_arms;ia++) {
        snr[ia] = 0.;
        if (spectro.lmin[ia]<373.8*(1+z) && 371.8*(1+z)<spectro.lmax[ia])
          snr[ia] = gsGetSNR_OII(&spectro,&obs,ia,z,1e-16,70.,0.0,1.0,spNoise[ia],snrType)
                    *sqrt((double)obs.n_exp)/1.6;
        snrtot += snr[ia]*snr[ia];
      }
      snrtot = sqrt(snrtot);   
      if (snrtot>snrmax16) snrmax16=snrtot;
    }
    mdlf = 1e-16/snrmax16*min_SNR*0.9; /* 0.9 is a safety factor */

    /* Modified by L.Dobos */
    /* Here we override previous r_eff values! */

    ngtot = 0;
    fp = fopen(InFileOII, "r");
    fq = fopen(OutFileOII, "w");
    while(fscanf(fp, "%ld %lg %lg %lg %lg %lg %lg",
      &id, &z, &obs.r_eff, &ROII, &FOII, &contOII, &sigma)!=EOF) {

      snrtot = 0.;
      if (FOII>=mdlf && obs.r_eff>=0.) for(ia=0;ia<spectro.N_arms;ia++) {
        snr[ia] = 0.;
        if (spectro.lmin[ia]<373.8*(1+z) && 371.8*(1+z)<spectro.lmax[ia])
          snr[ia] = gsGetSNR_OII(&spectro,&obs,ia,z,FOII,70.,contOII,ROII,spNoise[ia],snrType)
                    * sqrt((double)obs.n_exp);
        snrtot += snr[ia]*snr[ia];
      }
      snrtot = sqrt(snrtot);

      if (snrtot>=min_SNR) {
        fprintf(fq, "%7ld %5.3lf %5.3lf %11.5lE %9.4lf\n", id, z, obs.r_eff, FOII, snrtot);

        j = (long)floor((z-ZMIN_OII)/DZ_OII);
        if (j>=0 && j<NZ_OII) {ngal[j]++; ngtot++;}
      }
    }
    fclose(fp);
    fclose(fq);
    printf(" Done.\n");

    /* Report galaxy yields */
    printf("\n NUMBER OF AVAILABLE TARGETS: %7ld\n[i.e. objects where OII detected *if* targeted]\n\n", ngtot);
    printf(" zmin  zmax  ngal\n");
    for(j=0;j<NZ_OII;j++) {
      printf(" %5.3lf %5.3lf %7ld\n", ZMIN_OII+j*DZ_OII, ZMIN_OII+(j+1)*DZ_OII, ngal[j]);
    }
  }

  /* De-allocate noise vectors */
  gsFreeArmVectors(&spectro, spNoise);
  gsFreeArmVectors(&spectro, spSky);

  /* Deallocate magfiles */
  gsCopyMagfile(&inmag, &inmag2);

  return(0);
}