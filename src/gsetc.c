#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "libgsetc.h"

/* === MAIN PROGRAM === */

int main(void) {
  FILE *fp, *fq;
  int ia;
  long i, j, k;
  SPECTRO_ATTRIB spectro;
  OBS_ATTRIB obs;
  char FileName[256], OutFileNoise[256], OutFileSNR[256], OutFileSNRAB[256];
  char InFileOII[256], OutFileOII[256];
  int flag_reused;
  double lambda, z, t, fieldang, decent;
  double snr[MAXARM], snrtot, Aeff;
  double **spNoise;
  double **spSky;
  long id;
  double r_eff, ROII, FOII, contOII, sigma;
  double min_SNR = 0.;
  int n_exp, snrType;
  long ngal[NZ_OII], ngtot;
  double snrmax16, mdlf;
  double snrcont[MAXPIX];
  double snrcontcount[MAXPIX];
  double snrcontnoise[MAXPIX];
  double magcont[MAXPIX];
  double snctrans[MAXPIX];
  double samplefac[MAXPIX];
  int arm;
  long pix;
  double wav, innoise, inskymod;
  int proc, proc_tot;

  /* Added by Y.Moritani for input mag. file: 20150422*/
  char InFileMag[256], command[256], buf[256], dmy[32];
  FILE *in_pipe;
  double mag = 22.5;  /* default value of mag, used as a flag for inputfile*/

  /* Added by Y.Moritani for line SNR. file: 20150427*/
  char OutFileSNR2[256];
  double flux_emi, sigma_emi;

  /* Added by K.Yabe: 20160130 */
  double ref_input;

  double degrade;

  /* Tell us what flags are on */
  printf("Compiler flags:");
#ifdef DIFFRACTION_OFF
  printf(" -DDIFFRACTION_OFF");
#endif
#ifdef HGCDTE_SUTR
  printf(" -DHGCDTE_SUTR");
#endif
#ifdef NO_OBJ_CONTINUUM
  printf(" -DNO_OBJ_CONTINUUM");
#endif
#ifdef NATURAL_NLINES
  printf(" -DNATURAL_NLINES");
#endif
#ifdef MOONLIGHT_
  printf(" -DMOONLIGHT_");  
#endif
  printf("\n");

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

  gsReadSpectrographConfig(FileName, &spectro, degrade);

  /* Allocate noise vectors */
  spNoise=(double**)malloc((size_t)(spectro.N_arms*sizeof(double*)));
  for(ia=0; ia<spectro.N_arms; ia++) spNoise[ia] = (double*)malloc((size_t)(spectro.npix[ia]*sizeof(double)));
  spSky=(double**)malloc((size_t)(spectro.N_arms*sizeof(double*)));
  for(ia=0; ia<spectro.N_arms; ia++) spSky[ia] = (double*)malloc((size_t)(spectro.npix[ia]*sizeof(double)));

  /* Get observational conditions */
  /* Modified by Y.Moritnani -- 2016.02.16 */

  //printf("Enter observational conditions [hexadecimal code; suggested=11005]: ");
  if(scanf("%lx", &(obs.skytype))==EOF)
    obs.skytype = 0x10005;
/*
#if 0
  obs.skytype = 0x10005;
#endif
*/
  /* Input observing conditions */

  //printf("Enter seeing [arcsec FWHM @ lambda=800nm]: ");
  if(scanf("%lg", &(obs.seeing_fwhm_800))==EOF)
    obs.seeing_fwhm_800 = 0.8; 

  //printf("Enter zenith angle [degrees]: ");
  if(scanf("%lg", &(obs.zenithangle))==EOF)
    obs.zenithangle=45.;

  //printf("Enter Galactic dust column [magnitudes E(B-V)]: ");
  if(scanf("%lg", &(obs.EBV))==EOF)
    obs.EBV = 0.00;

  //printf("Enter field angle [degrees]: ");
  if(scanf("%lg", &fieldang)==EOF)
    fieldang = 0.675;

  //printf("Enter fiber astrometric offset [arcsec]: ");
  if(scanf("%lg", &decent)==EOF)
    decent = 0.00;

  /* Moon conditions -- set Moon below horizon unless otherwise indicated */
  obs.lunarZA = 135.;
  obs.lunarangle = 90.;
  obs.lunarphase = 0.25;
#ifdef MOONLIGHT_
  //printf("Enter Moonlight conditions:\n");
  //printf("  Angle (Moon-Zenith) [deg]: ");
  if(scanf("%lg", &(obs.lunarZA))==EOF)
    obs.lunarZA=30.;
  //printf("  Angle (Moon-target) [deg]: ");
  if(scanf("%lg", &(obs.lunarangle))==EOF)
    obs.lunarangle=60.;
  //printf("  Lunar phase [0=New, 0.25=quarter, 0.5=full]: ");
  if(scanf("%lg", &(obs.lunarphase))==EOF)
    obs.lunarphase=0.;
#endif
    
  /* Exposure time and systematics */
  //printf("Enter time per exposure [s]: ");
  if(scanf("%lg", &t)==EOF)
    t=450.;
  //printf("Enter number of exposures: ");
  if(scanf("%d", &n_exp)==EOF)
    n_exp=8;
  //printf("Enter systematic sky subtraction floor [rms per 1D pixel]: ");
  if(scanf("%lg", &(spectro.sysfrac))==EOF)
    spectro.sysfrac=0.01;
  spectro.sysfrac *= sqrt((double)n_exp); /* Prevent this from averaging down */
  //printf("Enter diffuse stray light [fraction of total]: ");
  if(scanf("%lg", &(spectro.diffuse_stray))==EOF)
    spectro.diffuse_stray=0.02;

  /* Modified by Y.Moritnani -- 2016.02.16 : end */

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
    sprintf(command,"wc %s", InFileMag);
    if ((in_pipe = popen(command,"r")) == NULL){
      exit(1);
    }
    if(fgets(buf,128, in_pipe)!=NULL)
      sscanf(buf, "%d %s", &num_inmag, dmy);

    lambda_inmag = (double*) malloc(sizeof(double) * num_inmag);
    mag_inmag = (double*) malloc(sizeof(double) * num_inmag);
    ia=0;
    fp=fopen(InFileMag,"r");
    while(fgets(buf, 256, fp) != NULL){
      if(strncmp(buf,"#",1) == 0) continue;
      sscanf(buf, "%lf %lf", lambda_inmag+ia, mag_inmag+ia);
      ia++;
    }
    num_inmag = ia;
    mag=-99.9; /* used as a flag for inputfile*/
  /* Modified by K.Yabe 20160703 */
    lambda_inmag2 = (double*) malloc(sizeof(double) * (num_inmag+2));
    mag_inmag2 = (double*) malloc(sizeof(double) * (num_inmag+2));
    if (lambda_inmag[0]>300.0){
      lambda_inmag2[0] = 300.0;
      mag_inmag2[0] = 99.9;
    }
    else {
      lambda_inmag2[0] = lambda_inmag[0]-1.0;
      mag_inmag2[0] = mag_inmag[0];
    }
    ia=0;
    while(ia<num_inmag){
      lambda_inmag2[ia+1] = lambda_inmag[ia];
      mag_inmag2[ia+1] = mag_inmag[ia];
      if (mag_inmag[ia]<=0.0){
        mag_inmag2[ia+1] = 99.9;
      }
      else {
        mag_inmag2[ia+1] = mag_inmag[ia];
      }
      ia++;
    }
    if (lambda_inmag[ia-1]<1300.0){
      lambda_inmag2[ia+1] = 1300.0;
      mag_inmag2[ia+1] = 99.9;
    }
    else{
      lambda_inmag2[ia+1] = lambda_inmag[ia-1]+1.0;
      mag_inmag2[ia+1] = mag_inmag[ia-1];
    }
    num_inmag2 = ia + 2;
/*    printf("%d %d\n", num_inmag, num_inmag2); */
    fclose(fp);
  }
  /* Added by Y.Moritani for input mag. file: 20150422 : end*/
  //printf("Enter the effective radius of the galaxy [arcsec]:\n");
  if(scanf("%lg", &ref_input)==EOF)
    ref_input=0.3;
  /* Added by K.Yabe */

  /* Encircled energy in fiber */
  printf("Fiber aperture factor [@800nm, r_eff=%.2lf\"(exp)] = %10.8lf\n", ref_input, gsGeometricThroughput(&spectro, &obs, 800, ref_input, decent, fieldang, 0x0));
  printf("Fiber aperture factor [@800nm,     point source] = %10.8lf\n", gsGeometricThroughput(&spectro, &obs, 800, 0.0, decent, fieldang, 0x0));

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
      gsGetNoise(&spectro,&obs,ia,fieldang,spNoise[ia],spSky[ia],t,0x0);
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
          snr[ia] = gsGetSNR_OII(&spectro,&obs,ia,z,flux_emi,sigma_emi,ref_input,0.,1.,decent,fieldang,spNoise[ia],t,0x0,snrType)
                    *sqrt((double)n_exp);
        snrtot += snr[ia]*snr[ia];
      }
      snrtot = sqrt(snrtot);

      Aeff = 0.;
      for(ia=0;ia<spectro.N_arms;ia++)
        if (spectro.lmin[ia]<372.845*(1+z) && 372.845*(1+z)<spectro.lmax[ia])
          Aeff += gsAeff(&spectro,&obs,ia,372.845*(1+z),fieldang);

      fprintf(fp, "%6.4lf %7.2lf %7.2lf %8.6lf %8.5lf",
        z, 372.71*(1+z), 372.98*(1+z), gsGeometricThroughput(&spectro, &obs, 372.845*(1+z), ref_input, decent, 0, 0x0), Aeff
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
          snr[ia] = gsGetSNR_Single(&spectro,&obs,ia,mag,345.5*(1+z),flux_emi,sigma_emi,ref_input,decent,fieldang,spNoise[ia],t,0x0,0)*sqrt((double)n_exp);
        snrtot += snr[ia]*snr[ia];
      }
      snrtot = sqrt(snrtot);

      Aeff = 0.;
      for(ia=0;ia<spectro.N_arms;ia++)
        if (spectro.lmin[ia]<345.5*(1+z) && 345.5*(1+z)<spectro.lmax[ia])
          Aeff += gsAeff(&spectro,&obs,ia,345.5*(1+z),fieldang);

/*      fprintf(fp, "%6.4lf %7.2lf %8.6lf %8.5lf", z, 345.5*(1+z), gsGeometricThroughput(&spectro, &obs, 345.5*(1+z), ref_input, decent, 0, 0x0), Aeff); */
      fprintf(fp, "%7.2lf %8.6lf %8.5lf", 345.5*(1+z), gsGeometricThroughput(&spectro, &obs, 345.5*(1+z), ref_input, decent, fieldang, 0x0), Aeff);
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
      /* Modified by Y. Moritani for input mag. file: 20150422 :*/
      /* Modified by K. Yabe for counts output: 20150525 :*/
      //gsGetSNR_Continuum(&spectro,&obs,ia,22.5,0.0,decent,fieldang,spNoise[ia],t,0x0,snrcont);
      gsGetSNR_Continuum(&spectro,&obs,ia,mag,ref_input,decent,fieldang,spNoise[ia],t,0x0,snrcont,snrcontcount,snrcontnoise,magcont,snctrans,samplefac);
      for(j=0;j<spectro.npix[ia];j++) {
        fprintf(fp, "%1d %4ld %9.3lf %8.4lf %11.5lE %11.5lE %11.5lE %11.5lE %11.5lE %11.5lE  %11.5lE\n",
		spectro_arm(&spectro, ia), j, spectro.lmin[ia]+spectro.dl[ia]*j,snrcont[j]*sqrt((double)n_exp),snrcontcount[j],spNoise[ia][j],snrcontnoise[j],magcont[j],snctrans[j],samplefac[j],spSky[ia][j]);
      }
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
          snr[ia] = gsGetSNR_OII(&spectro,&obs,ia,z,1e-16,70.,0.0,0.0,1.0,decent,fieldang,spNoise[ia],t,0x0,snrType)
                    *sqrt((double)n_exp)/1.6;
        snrtot += snr[ia]*snr[ia];
      }
      snrtot = sqrt(snrtot);   
      if (snrtot>snrmax16) snrmax16=snrtot;
    }
    mdlf = 1e-16/snrmax16*min_SNR*0.9; /* 0.9 is a safety factor */

    ngtot = 0;
    fp = fopen(InFileOII, "r");
    fq = fopen(OutFileOII, "w");
    while(fscanf(fp, "%ld %lg %lg %lg %lg %lg %lg",
      &id, &z, &r_eff, &ROII, &FOII, &contOII, &sigma)!=EOF) {

      snrtot = 0.;
      if (FOII>=mdlf && r_eff>=0.) for(ia=0;ia<spectro.N_arms;ia++) {
        snr[ia] = 0.;
        if (spectro.lmin[ia]<373.8*(1+z) && 371.8*(1+z)<spectro.lmax[ia])
          snr[ia] = gsGetSNR_OII(&spectro,&obs,ia,z,FOII,70.,r_eff,contOII,ROII,decent,fieldang,spNoise[ia],t,0x0,snrType)
                    * sqrt((double)n_exp);
        snrtot += snr[ia]*snr[ia];
      }
      snrtot = sqrt(snrtot);

      if (snrtot>=min_SNR) {
        fprintf(fq, "%7ld %5.3lf %5.3lf %11.5lE %9.4lf\n", id, z, r_eff, FOII, snrtot);

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
  for(ia=0; ia<spectro.N_arms; ia++) free((char*)(spNoise[ia]));
  free((char*)spNoise);
  for(ia=0; ia<spectro.N_arms; ia++) free((char*)(spSky[ia]));
  free((char*)spSky);
  return(0);
}