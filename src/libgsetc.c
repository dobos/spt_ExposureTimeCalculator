#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include "libgsetc.h"

/* Sky model data */
#include "modeldata.h"

/* --- SPECIAL FUNCTIONS --- */
  
/* Bessel functions, J0 and J1:
 * These are needed to convert between real and Fourier space.
 * The |x|>3 expansions are formulae 9.4.3,9.4.6 of Abramowitz &
 * Stegun, good to several parts in 10^8 (i.e. good enough for
 * PSF prediction work). At |x|<3 the integral definition
 * is used.
 */
  
double getJ0(double x) {
  double f,t,u;
  int i;

  x = fabs(x);
    
  /* Large value - Abramowitz & Stegun */
  if (x>3) {
    u=3./x;
    f = 0.79788456 - 0.00000077*u - 0.00552740*u*u - 0.00009512*u*u*u
       + u*u*u*u*(0.00137234-0.00072805*u+0.00014476*u*u);
    t = x - 0.78539816 - 0.04166397*u - 0.00003954*u*u + 0.00262573*u*u*u
       + u*u*u*u*(-0.00054125 -0.00029333*u +0.00013558*u*u);
    return(f/sqrt(x)*cos(t));
  }
  
  /* Small value - Abramowitz & Stegun */
  f = 0.;
  for(i=0;i<50;i++) f+=cos(x*cos((i+0.5)*M_PI/50.));
  return(0.02*f);   
}
   
double getJ1(double x) {
  double f,t,u,s;
  int i;
  
  /* Flip sign for negative values */
  s=1; if (x<0) {x=-x; s=-1;}
   
  /* Large value - Abramowitz & Stegun */
  if (x>3) {
    u=3./x;
    f = 0.79788456 + 0.00000156*u + 0.01659667*u*u + 0.00017105*u*u*u
       + u*u*u*u*(-0.00249511 + 0.00113653*u - 0.00020033*u*u);
    t = x - 2.35619449 + 0.12499612*u + 0.0000565*u*u - 0.00637879*u*u*u
       + u*u*u*u*(0.00074348 + 0.00079824*u - 0.00029166*u*u);
    return(s*f/sqrt(x)*cos(t));
  }

  /* Small value - Abramowitz & Stegun */
  f = 0.;
  for(i=0;i<50;i++) f+=cos(x*sin((i+0.5)*M_PI/50.)-((i+0.5)*M_PI/50.));
  return(0.02*s*f);
}
 
/* Error function */
double geterf(double x) {
  
  int j;
  double s, term, erfabsx, erfcx, y, dy, u;
  
  s=1;
  if (x<0) {  
    s=-1; x=-x;
  }

  /* For values greater than 6, the result for erf will be unity to machine precision */
  if (x>6) return(s);
       
  /* Taylor expansion for smallest values */
  if (x<1.5) {
    erfabsx = 0.;
    term = 2./sqrt(M_PI)*x;
    for(j=0;j<=24;j++) {
      erfabsx += term/(2*j+1);
      term *= -x*x/(j+1);
    }
    return(s*erfabsx);
  }
   
  /* Compute erfc(x) by transforming the complementary error function integral:
   * erfc(x) = 2/sqrtpi * int_x^infty e^(-z^2) dz
   * transform: z = x + exp(y)/x, dz = exp(y)/x * dy
   * erfc(x) = 2/(sqrtpi*x) * int_-infty^infty e^{y-[x+exp(y)/x]^2} dx
   * The integrand rises rapidly, peaks near y~-0.7, and then falls.
   * Convergence is exponential or faster in both directions.
   * Current integration range -30<y<3 chosen for accuracy at 1.5<x<6.
   */
  erfcx = 0.;
  dy = 0.01;
  for(j=-3000;j<300;j++) {
    y=j*dy;
    u = x+exp(y)/x;
    erfcx += exp(y-u*u);
  }
  erfcx *= 2./sqrt(M_PI)/x*dy;
  erfabsx=1-erfcx;
  
  return(s*erfabsx);
}

/* --- CONVERSION FUNCTIONS --- */

/* Index of refraction of air -- Edlen 1953 formula */
double gs_n_air(double lambda_vac) {
  return(1.000064328 + 0.0294981/(146-1e6/lambda_vac/lambda_vac) + 0.0002554/(41-1e6/lambda_vac/lambda_vac));
}

/* Vacuum -> air conversion */
double gs_vac2air(double lambda_vac) {
  return(lambda_vac/gs_n_air(lambda_vac));
}

/* Air -> vacuum conversion */
double gs_air2vac(double lambda_air) {
  double lambda_vac;
  int i;

  lambda_vac = lambda_air;
  for(i=0;i<10;i++) {
    lambda_vac += lambda_air - gs_vac2air(lambda_vac);
  }
  return(lambda_vac);
}

/* --- MATERIAL PROPERTIES --- */

/* Silicon absorption as a function of temperature and wavelength.
 * lambda in nm, T in Kelvin, output in microns
 *
 * From: Rajkanan, Singh, and Shewchun, Solid-State Electronics 22, 793-795, 1979.
 *
 * Loosely based on the IDL code by D Groom
 */
double gsOP_Si_abslength(double lambda, double T) {

  int i,j;
  double hnu, E_gd, E_g[2], alpha, de;

  double beta = 7.021e-4;              /* eV/K        */
  double k = 8.617e-5;                 /* eV/K        */
  double gamma = 1108.;                /* K           */
  double E_g0[] = {1.1557, 2.5};       /* eV          */
  double E_gd0 = 3.2;                  /* eV          */
  double E_p[] = {1.827e-2, 5.773e-2}; /* eV          */
  double C[] = {5.5, 4.0};             /* 1           */
  double A[] = {323.1, 7237.};         /* cm^-1 eV^-2 */
  double A_d = 1.052e6;                /* cm^-1 eV^-2 */

  lambda /= 1e9; /* Convert to meters */

  /* Legal range: 2e-7 .. 1.1e-6 meters */
  if (lambda<2e-7 || lambda>1.1e-6) {
    gsError("Error: Si data out of range.");
  }

  hnu = 1.239842e-6/lambda; /* in eV */
  E_gd = E_gd0 - beta*T*T/(T+gamma);
  for(i=0;i<2;i++) E_g[i] = E_g0[i] - beta*T*T/(T+gamma);
  alpha = 0.;
  if (hnu>E_gd) alpha += A_d*sqrt(hnu-E_gd);
  for(i=0;i<2;i++) for(j=0;j<2;j++) {
    de = hnu - E_g[j] + E_p[i];
    if (de>0) alpha += C[i]*A[j]*de*de/( exp(E_p[i]/k/T)-1 );
    de = hnu - E_g[j] - E_p[i];
    if (de>0) alpha += C[i]*A[j]*de*de/( 1-exp(-E_p[i]/k/T) );
  }

  return(1e4/alpha);
}

/* Index of refraction of silicon as a function of wavelength (in nm).
 *
 * (Right now the temperature is vestigial: assumes room temperature. Probably
 * OK for real part.)
 *
 * Interpolation based on: http://snap.lbl.gov/ccdweb/ccd_data/complex_index.dat
 *
 * 1 Interpolated table of optical constants for silicon. Don Groom 10 Dec 98
 * 2 See SPIE 1999 reference for details. Complex n_c = n + i k.
 * 3 For lambda < 750 nm data is from D. F. Edwards in the Handbook of Optical
 * 4 Constants of Solids (1985). Reflectivity and absorption curves from
 * 5 Janesick's SPIE course notes have been used for redward extention.
 * 6 Small differences between tabulated k and \lambda/4\pi\ell are due
 * 7 to curve smoothing. Significance on n and k overstated for interpolation.
 * 8 f='(f10.1,f10.4,f10.4,f10.4,e12.3)'
 */
double gsOP_Si_indexreal(double lambda, double T) {

  double si_table[] = {
    0.9317, 1.0031, 1.0811, 1.1490, 1.2253, 1.3295, 1.4848, 1.5861, 1.5833, 1.5690, 1.5800, 1.6170, 1.6834, 1.7971, 2.0176, 2.4043, 2.9198, 3.5317, 4.3658, 4.8761, 
    5.0036, 5.0203, 5.0100, 5.0096, 5.0234, 5.0551, 5.0983, 5.1548, 5.2220, 5.3086, 5.4358, 5.6569, 6.0519, 6.5467, 6.8110, 6.7364, 6.4693, 6.1850, 5.9438, 5.7404, 
    5.5700, 5.4254, 5.2960, 5.1868, 5.0890, 5.0024, 4.9235, 4.8522, 4.7873, 4.7289, 4.6738, 4.6230, 4.5759, 4.5317, 4.4918, 4.4543, 4.4202, 4.3866, 4.3553, 4.3253, 
    4.2975, 4.2716, 4.2460, 4.2224, 4.2000, 4.1787, 4.1586, 4.1378, 4.1197, 4.1017, 4.0844, 4.0682, 4.0526, 4.0376, 4.0226, 4.0092, 3.9954, 3.9825, 3.9700, 3.9585, 
    3.9473, 3.9367, 3.9262, 3.9156, 3.9058, 3.8955, 3.8865, 3.8776, 3.8684, 3.8594, 3.8512, 3.8435, 3.8362, 3.8285, 3.8208, 3.8134, 3.8066, 3.8005, 3.7946, 3.7888, 
    3.7831, 3.7774, 3.7712, 3.7660, 3.7617, 3.7566, 3.7514, 3.7474, 3.7430, 3.7379, 3.7333, 3.7289, 3.7250, 3.7212, 3.7176, 3.7139, 3.7093, 3.7048, 3.7008, 3.6968, 
    3.6925, 3.6881, 3.6848, 3.6815, 3.6778, 3.6742, 3.6719, 3.6703, 3.6686, 3.6670, 3.6654, 3.6637, 3.6621, 3.6605, 3.6589, 3.6572, 3.6556, 3.6540, 3.6523, 3.6503, 
    3.6482, 3.6460, 3.6439, 3.6418, 3.6397, 3.6375, 3.6354, 3.6333, 3.6311, 3.6290, 3.6269, 3.6248, 3.6233, 3.6225, 3.6217, 3.6203, 3.6183, 3.6163, 3.6144, 3.6125, 
    3.6106, 3.6091, 3.6076, 3.6061, 3.6046, 3.6031, 3.6015, 3.6000, 3.5984, 3.5968, 3.5953, 3.5937, 3.5922, 3.5906, 3.5890, 3.5875, 3.5859, 3.5844, 3.5828, 3.5812, 
    3.5797};
  double x, xfrac;
  int xint;

  lambda /= 1e9; /* Convert to meters */

  /* Legal range: 2e-7 .. 1.1e-6 meters */
  if (lambda<2e-7 || lambda>1.1e-6) {
    gsError("Error: Si data out of range.");
  }

  /* Get position in table to read */
  x = (lambda-2e-7)/5e-9;
  xint = (int)floor(x);
  if (xint<=  0) xint=  1;
  if (xint>=179) xint=178;
  xfrac = x - xint;

  return(
    -xfrac*(xfrac-1)*(xfrac-2)/6.*si_table[xint-1]
    +(xfrac*xfrac-1.)*(xfrac-2)/2.*si_table[xint]
    -xfrac*(xfrac+1)*(xfrac-2)/2.*si_table[xint+1]
    +xfrac*(xfrac*xfrac-1)/6.*si_table[xint+2]
  );
}

/* --- ROUTINES TO COMPUTE THE FOREGROUND ABSORPTION --- */
    
/* Ratio of extinction A_lambda to reddening E(B-V) for Milky Way dust.
 * Input lambda in nm.
 *
 * Valid range: 0.1--10.0 um
 *
 * Uses lookup table based on:
 * ftp://ftp.astro.princeton.edu/draine/dust/mix/kext_albedo_WD_MW_3.1_60_D03.all
 * Retrieved 03/18/2011
 */
double gsGalactic_Alambda__EBV(double lambda) {
  
  int il;
  double lset, lf;
  double norm[] = {
    0.24174,  0.25504,  0.26708,  0.26392,  0.23976,  0.20540,  0.17650,  0.15484,  0.13199,  0.10599,
    0.08604,  0.07340,  0.06302,  0.05489,  0.05107,  0.05005,  0.05041,  0.05151,  0.05313,  0.05554,
    0.06009,  0.06392,  0.06130,  0.06117,  0.06239,  0.06417,  0.06629,  0.06860,  0.07123,  0.07400,
    0.07702,  0.08018,  0.08361,  0.08723,  0.09105,  0.09506,  0.09934,  0.10388,  0.10862,  0.11356,
    0.11876,  0.12416,  0.12989,  0.13581,  0.14207,  0.14858,  0.15543,  0.16261,  0.17156,  0.17781,
    0.18571,  0.19401,  0.20276,  0.21192,  0.22140,  0.23127,  0.24147,  0.25214,  0.26313,  0.27459,
    0.28650,  0.29875,  0.31152,  0.32469,  0.33831,  0.35240,  0.36708,  0.38216,  0.39776,  0.41389,
    0.43055,  0.44779,  0.46557,  0.48394,  0.50296,  0.52265,  0.54292,  0.56386,  0.58552,  0.60790,
    0.63094,  0.65477,  0.67939,  0.70441,  0.73074,  0.75839,  0.78670,  0.81567,  0.84529,  0.87689,
    0.90915,  0.94273,  0.97762,  1.01382,  1.05201,  1.09151,  1.13232,  1.17512,  1.21922,  1.26531,
    1.31336,  1.36340,  1.41475,  1.46741,  1.52205,  1.57867,  1.63594,  1.69585,  1.75708,  1.82028,
    1.88545,  1.95260,  2.02107,  2.09217,  2.16524,  2.24029,  2.31731,  2.39631,  2.47795,  2.56155,
    2.64714,  2.73469,  2.82423,  2.91573,  3.00987,  3.10533,  3.20342,  3.30349,  3.40487,  3.50823,
    3.61290,  3.71955,  3.82752,  3.93680,  4.04608,  4.15668,  4.26728,  4.37854,  4.49111,  4.60369,
    4.71692,  4.83081,  4.94470,  5.05991,  5.17643,  5.29361,  5.41145,  5.53061,  5.65109,  5.77354,
    5.89862,  6.02567,  6.15339,  6.26596,  6.38907,  6.52732,  6.68203,  6.87294,  7.09019,  7.33377,
    7.64319,  7.96577,  8.45293,  8.95326,  9.47334,  9.94075, 10.22383, 10.19750,  9.92100,  9.52600,
    9.11126,  8.75576,  8.47268,  8.26201,  8.11718,  8.01843,  7.96577,  7.95260,  7.96577,  8.02502,
    8.11718,  8.24226,  8.39368,  8.55826,  8.74259,  8.92034,  9.10467,  9.31534,  9.54575,  9.82883,
   10.11192, 10.41475, 10.72416, 11.05332, 11.40224, 11.77749, 12.17248, 12.61356, 13.06781, 13.52205,
   13.92363};

  lambda /= 1e3; /* convert nm --> um */

  if (lambda<.1 || lambda>10) {  
    gsError("Error: Wavelength lambda = %12.5lE microns out of range for dust extinction law.", lambda);
  }
   
  /* Find place to linearly interpolate */
  lset = 100./log(10.) * log(10./lambda);
  il = (int)floor(lset);
  if (il<0) il=0;
  if (il>199) il=199;
  lf = lset-il;
    
  return(norm[il] + lf*(norm[il+1]-norm[il]));
}

/* Continuum atmospheric opacity in magnitudes per airmass, as a function
 * of wavelength lambda (in nm)
 */
double gsAtmContOp(OBS_ATTRIB *obs, double lambda) {
  double k;

  k = 0.113; /* V band opacity -- placeholder !!! */

  /* Opacity model */
  switch((obs->skytype>>8) & 0xf) {

    case 0x0:
      /* Median extinction on Mauna Kea -- Gemini website
       * Boulade 1988, CFHT Bulletin, 19, 16, for l<400nm; CFHT Observers Manual for l>400nm
       */
      k=0.05;
      if (lambda<900) k=0.07-0.02*(lambda-800.)/100.;
      if (lambda<800) k=0.10-0.03*(lambda-700.)/100.;
      if (lambda<700) k=0.11-0.01*(lambda-650.)/ 50.;
      if (lambda<650) k=0.11-0.00*(lambda-600.)/ 50.;
      if (lambda<600) k=0.12-0.01*(lambda-550.)/ 50.;
      if (lambda<550) k=0.13-0.01*(lambda-500.)/ 50.;
      if (lambda<500) k=0.17-0.04*(lambda-450.)/ 50.;
      if (lambda<450) k=0.25-0.08*(lambda-400.)/ 50.;
      if (lambda<400) k=0.30-0.05*(lambda-380.)/ 20.;
      if (lambda<380) k=0.37-0.07*(lambda-360.)/ 20.;
      if (lambda<360) k=0.51-0.14*(lambda-340.)/ 20.;
      if (lambda<340) k=0.82-0.31*(lambda-320.)/ 20.;
      if (lambda<320) k=1.37-0.55*(lambda-310.)/ 10.;
      if (lambda<310) k=1.37;
      break;

    default:
      gsError("Error: gsAtmContOp: illegal opacity model.");
      break;
  }

  return(k);
}


/* Atmospheric transmission, as a function of wavelength lambda (in nm)
 * and observing conditions.
 */
double gsAtmTrans(OBS_ATTRIB *obs, double lambda) {
  double k, trans;
  double x, xfrac;
  long xint;

  k = gsAtmContOp(obs,lambda);
  trans = pow(10., -0.4*k/cos(obs->zenithangle*DEGREE));

  /* Absorption lines */
  switch((obs->skytype>>12) & 0xf) {
    
    case 0x0:
      /* Currently from Kitt Peak transmission spectrum */
      x = (lambda-500)/0.025;
      xint = (long)floor(x);
      xfrac = x-xint;
      if (xint>=0 && xint<=39998) {
        /* Optical */
        trans *= AtmTransKP[xint] * (1-xfrac) + AtmTransKP[xint+1] * xfrac;
      }
      if (xint>39998) {
        /* Infrared */
        gsError("Error: transmission model not valid in the IR.");
      }
      break;

    case 0x1:
      /* Kitt Peak at lambda<900nm; for lambda>900nm switch to Gemini simulated spectrum
       * at 3mm precipitable water.
       */
      if (lambda>900) {
        x = (lambda-900)/0.02;
        xint = (long)floor(x);
        xfrac = x-xint;
        if (xint<0) {xint=0; xfrac=0.;}
        if (xint<30000) {
          /* Far-red -- 3 mm water column */
          trans *= MKTrans_3mm[xint] * (1-xfrac) + MKTrans_3mm[xint+1] * xfrac;
        } else {
          /* Infrared, beyond 1500nm */
          gsError("Error: transmission model not valid in the IR.");
        }
      } else {
        x = (lambda-500)/0.025;
        xint = (long)floor(x);
        xfrac = x-xint;
        if (xint>=0)
          trans *= AtmTransKP[xint] * (1-xfrac) + AtmTransKP[xint+1] * xfrac;
      }
      break;

    default:
      gsError("Error: Unrecognized atmospheric line absorption model.");
      break;
  }

  return(trans);
}

/* --- ROUTINES TO CALCULATE THROUGHPUT & RELATED QUANTITIES --- */

/* Computes the geometric throughput, or encircled energy in a fiber.
 * Requires attributes and:
 *   lambda (wavelength in nm)
 *   r_eff (target effective radius in arcsec; 0 for star)
 *   decent (decenter of fiber from target in arcsec)
 *   fieldang (field angle in degrees; 0 on axis)
 *
 * Returns 0 if failed.
 */
double gsGeometricThroughput(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, double lambda) {

  int i;
  long iu, Nu;
  double u, k, du, ftilde, Gtilde, EE;
  double R, sigma, uscale, rs, theta_D, EFL;

  /* The parameters */
  Nu = 50;
  du = 0.024;

  /* Smearing in arcsec, and the EFL */
  i = floor(4*obs->fieldangle/spectro->rfov);
  if (i<0) {
    sigma = spectro->rms_spot[0];
    EFL = spectro->EFL[0];
  } else {
    if (i>=MAXEFL-1) {
      sigma = spectro->rms_spot[MAXEFL-1];
      EFL = spectro->EFL[MAXEFL-1];
    } else {
      sigma = spectro->rms_spot[i]
              + (spectro->rms_spot[i+1]-spectro->rms_spot[i])
              * (4*obs->fieldangle/spectro->rfov-i);
      EFL = spectro->EFL[i]
              + (spectro->EFL[i+1]-spectro->EFL[i])
              * (4*obs->fieldangle/spectro->rfov-i);
    }
  }
  sigma *= ARCSEC_PER_URAD / EFL; 

  /* Fiber radius in arcsec */
  R = spectro->fiber_ent_rad / EFL * ARCSEC_PER_URAD;

  /* Seeing MTF = 1/e @ u=uscale */
  uscale = 0.465/obs->seeing_fwhm_800*pow(lambda/800.,0.2);

  /* Galaxy scale length */
  rs = obs->r_eff/RAT_HL_SL_EXP;

  /* The integral */
  EE = 0.;
  for(iu=0;iu<Nu;iu++) {
    u = (iu+0.5)*du;
    k = 2.*M_PI*u;

    /* Telescope PSF.
     * Treats diffraction using only the leading term in u, which should be okay
     * for telescopes where seeing dominates and the only importance of the diffraction
     * is to scatter a small fraction of light into the far wings.
     */
    theta_D = 0.001*lambda/spectro->D_outer/(1.-spectro->centobs) * ARCSEC_PER_URAD;
#ifdef DIFFRACTION_OFF
    theta_D = 0.;
#endif
    Gtilde = exp(-k*k*sigma*sigma/2.-pow(u/uscale,1.666666666666666666666666667))
             * getJ0(k*obs->decent)
             * exp(-4./M_PI*u*theta_D);

    /* Galaxy profile */
    ftilde = pow(1. + k*rs*k*rs, -1.5);

    EE += R*2.*M_PI*du*getJ1(2.*M_PI*u*R)*ftilde*Gtilde;
  }
  return(EE);
}

/* Computes the effective area of the system including all
 * throughput factors from the telescope through detector (but
 * not the atmosphere).
 */
double gsAeff(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm) {
  double Aeff, Vig;
  int i;

  /* Geometric area */
  Aeff = M_PI/4. * spectro->D_outer * spectro->D_outer * (1-spectro->centobs*spectro->centobs);

  /* Vignetting */
  i = floor(4*obs->fieldangle/spectro->rfov);  
  if (i<0) {
    Vig = spectro->vignette[0];
  } else {
    if (i>=4) { 
      Vig = spectro->vignette[4];
    } else {
      Vig = spectro->vignette[i]
              + (spectro->vignette[i+1]-spectro->vignette[i])
              * (4*obs->fieldangle/spectro->rfov-i);
    }
  }
  Aeff *= Vig;

  return Aeff;
}

/* Interpolate the troughtput (excluding atmosphere & fiber geom)
 */
double gsThroughput(SPECTRO_ATTRIB *spectro, int i_arm, double lambda) {
  double Thr, fr;
  int imin, imax;
  int ti;

  /* Throughput */
  imin = spectro->istart[i_arm];
  imax = spectro->istart[i_arm+1];
  if (lambda<=spectro->l[imin]) {
    Thr = spectro->T[imin];
  } else {
    if (lambda>=spectro->l[imax-1]) {
      Thr = spectro->T[imax-1];
    } else {
      ti=imin;
      while (lambda>spectro->l[ti+1]) ti++;
      fr = (lambda-spectro->l[ti])/(spectro->l[ti+1]-spectro->l[ti]);
      Thr = spectro->T[ti] + (spectro->T[ti+1]-spectro->T[ti])*fr;
    }
  }
  return Thr;
}

/* Computes the 1D Fourier transform (in the dispersion direction) of the spectrograph PSF
 * at wavenumber u (units: cycles/pixel). We treat only the real part (i.e. ignore skewed
 * PSFs).
 */
double gsSpectroMTF(SPECTRO_ATTRIB *spectro, int i_arm, double lambda, double u) {

  double D_spot, sigma;
  double mtf = 1.;
  int i, N;
  double depth, ddepth, contrib, numer, denom, mfp, nSi, rspot, d0, arg;

  /* Fiber size */
  D_spot = spectro->diam[i_arm]/spectro->pix[i_arm];
  mtf *= fabs(D_spot*u)>1e-6? 2.*getJ1(M_PI*D_spot*u)/(M_PI*D_spot*u): 1.;

  /* Pixelization */
  mtf *= fabs(u)>1e-9? sin(M_PI*u)/(M_PI*u): 1.;

  /* Spot size */
  sigma = spectro->rms_cam[i_arm]/spectro->pix[i_arm];
  mtf *= exp(-2.*M_PI*M_PI*sigma*sigma*u*u);

  /* Defocus in the detector -- Si only */
  if (spectro->Dtype[i_arm]==0) {
    N = 3;
    ddepth = spectro->thick[i_arm]/N;
    numer = denom = 0.;
    mfp = gsOP_Si_abslength(lambda, spectro->temperature[i_arm]);
    nSi = gsOP_Si_indexreal(lambda, spectro->temperature[i_arm]);
    d0 = mfp<1e4*ddepth? mfp * ( 1 - (1+ddepth/mfp)*exp(-ddepth/mfp) )/( 1 - exp(-ddepth/mfp) ): 0.5*ddepth;
    for(i=0;i<2*N;i++) {
      depth = d0 + ddepth*i;
      contrib = exp(-ddepth*i/mfp) * (depth>spectro->thick[i_arm]?0.3:1.0);
      rspot = depth/nSi/2./spectro->fratio[i_arm]/spectro->pix[i_arm];
      arg = 2*M_PI*rspot*u;
      denom += contrib;

      /* MTF modeled using approximation to integral,
       * MTF[single depth] = <sin^2 phi cos(arg * cos phi)> averaged over phi
       * Use steps of pi/6 here.
       */
      numer += contrib * (0.1666666667*cos(0.866025404*arg) + 0.5*cos(0.5*arg) + 0.3333333333);
    }
    mtf *= numer/denom;
  }

  /* Scattering from grating */
  mtf *= exp(-lambda/spectro->dl[i_arm]/spectro->nline[i_arm]*fabs(u));

  return(mtf);
}

/* For a feature centered at position pos (in pixels), computes the fraction of radiation
 * in each l-pixel (from 0 .. N-1); returns to fr[0..N-1]. Adds additional smearing of
 * Gaussian width sigma (in pixels).
 */
void gsSpectroDist(SPECTRO_ATTRIB *spectro, int i_arm, double lambda,
  double pos, double sigma, int N, double *fr) {

  int ip;
  long iu, Nu;
  double u, du, mtf1d;

  for(ip=0;ip<N;ip++) fr[ip] = 0.;
  Nu = 1000; du = 0.005;
  for(iu=0;iu<Nu;iu++) {
    u = du*(iu+0.5);
    mtf1d = gsSpectroMTF(spectro,i_arm,lambda,u) * exp(-2.*M_PI*M_PI*sigma*sigma*u*u);
    for(ip=0;ip<N;ip++)
      fr[ip] += 2.*du*cos(2.*M_PI*u*(pos-ip))*mtf1d;
  }

  return;
}

/* Routine to compute the fraction of radiation landing in a bin of spectro->width pixels.
 * Treats only the radiation from the correct trace (tr=0) or from adjacent traces as well
 * (tr=1).
 */
double gsFracTrace(SPECTRO_ATTRIB *spectro, int i_arm, double lambda, int tr) {

  double *FR, sum;
  int i, j, N;

  N = spectro->width[i_arm];
  FR = (double*)malloc((size_t)(N*sizeof(double)));
  sum=0;
  for(j=-tr;j<=tr;j++) {
    gsSpectroDist(spectro,i_arm,lambda,0.5*(N-1)+j*spectro->sep[i_arm]/spectro->pix[i_arm],0,N,FR);
    for(i=0;i<N;i++) sum += FR[i];
  }
  free((char*)FR);
  return(sum);
}

/* --- ROUTINES TO COMPUTE THE SIGNAL AND THE NOISE --- */

double gsGetAirmass(OBS_ATTRIB* obs) {
  double airmass;

  /* Airmass computation */
  airmass = 1./sqrt(1.-0.96*sin(obs->zenithangle*DEGREE)*sin(obs->zenithangle*DEGREE));
  return airmass;
}

double gsGetFiberRadius(SPECTRO_ATTRIB* spectro, OBS_ATTRIB* obs, int i_arm) {
  int i;
  double EFL, rad;

  /* Fiber radius in arcsec */
  i = floor(4*obs->fieldangle/spectro->rfov);
  if (i<0) {
    EFL = spectro->EFL[0]; 
  } else {
    if (i>=MAXEFL-1) {
      EFL = spectro->EFL[MAXEFL-1];
    } else {
      EFL = spectro->EFL[i]
              + (spectro->EFL[i+1]-spectro->EFL[i])
              * (4*obs->fieldangle/spectro->rfov-i);
    }
  }
  rad = spectro->fiber_ent_rad/EFL * ARCSEC_PER_URAD;
  return rad;
}

double gsGetSampleFactor(SPECTRO_ATTRIB* spectro, int i_arm) {
  double sample_factor;

  /* Pre-factor for sampling the Poisson distribution */
  sample_factor = 1.0;
#ifdef HGCDTE_SUTR
  if (spectro->Dtype[i_arm]==1) sample_factor = 1.2;
#endif

  return sample_factor;
}

void gsAddSkyLines_UVES(SPECTRO_ATTRIB* spectro, OBS_ATTRIB* obs, int i_arm, double* Noise) {
  long iline, iref, j;
  double airmass, rad, sample_factor;
  double lambda, pos;
  double count;
  double FR[5*SP_PSF_LEN];
  
  airmass = gsGetAirmass(obs);
  rad = gsGetFiberRadius(spectro, obs, i_arm);
  sample_factor = gsGetSampleFactor(spectro, i_arm);

  for(iline=0;iline<gsSKY_UVES_NLINES;iline++) {
    lambda = gs_air2vac(gsSKY_UVES_LAMBDA[iline]);
    pos = (lambda-spectro->lmin[i_arm])/spectro->dl[i_arm];
    if (pos>-(SP_PSF_LEN/2) && pos<spectro->npix[i_arm]+SP_PSF_LEN/2-1) {
      /* This line is within the range of this spectrograph arm.
        * Need to get how the counts are distributed and then add this line.
        */

      /* Inputs are in units of 1e-12 erg/m^2/s/arcsec^2 --> need to do
        * appropriate conversion to counts in detector.
        */
      count = gsSKY_UVES_INT[iline] * lambda * 1e-12 * PHOTONS_PER_ERG_1NM
              * gsFracTrace(spectro,i_arm,lambda,1)
              * gsAeff(spectro,obs,i_arm) 
              * gsThroughput(spectro,i_arm,lambda)
              * obs->t_exp * M_PI * rad * rad;
      if (count<0) count=0;

      /* Rescale line counts by the airmass; UVES referenced to 1.1. Also by atmospheric
        * extinction curve.
        */
      count *= airmass/1.1  * exp(-gsAtmContOp(obs,lambda)*airmass/1.086);

      iref = (long)floor(pos-7.5);
      if (iref<0) iref=0;
      if (iref>spectro->npix[i_arm]-16) iref=spectro->npix[i_arm]-16;
      gsSpectroDist(spectro,i_arm,lambda,pos-iref,0,16,FR);
      for(j=0;j<16;j++)
        Noise[iref+j] += count*FR[j]*sample_factor;
    }
  }
}

void gsAddSkyLines_NIR(SPECTRO_ATTRIB* spectro, OBS_ATTRIB* obs, int i_arm, double* Noise) {
  long iline, iref, j;
  double airmass, rad, sample_factor;
  double lambda, pos;
  double count;
  double FR[5*SP_PSF_LEN];

  airmass = gsGetAirmass(obs);
  rad = gsGetFiberRadius(spectro, obs, i_arm);
  sample_factor = gsGetSampleFactor(spectro, i_arm);

  /* Now the NIR lines */
  for(iline=0;iline<N_IR_OH_LINE;iline++) {

    lambda = OHDATA[2*iline];
  
    pos = (lambda-spectro->lmin[i_arm])/spectro->dl[i_arm];    
    if (pos>-(SP_PSF_LEN/2) && pos<spectro->npix[i_arm]+SP_PSF_LEN/2-1) { 
      /* Inputs are in units of 1e-12 erg/m^2/s/arcsec^2 --> need to do
        * appropriate conversion to counts in detector.
        */
      count = OHDATA[2*iline+1] * lambda * 1e-12 * PHOTONS_PER_ERG_1NM
              * gsFracTrace(spectro,i_arm,lambda,1)
              * gsAeff(spectro,obs,i_arm)
              * gsThroughput(spectro,i_arm,lambda)
              * obs->t_exp * M_PI * rad * rad;
      if (count<0) count=0;
  
      /* Rescale line counts by the airmass; referenced to 1.0. Also by atmospheric
        * extinction curve, and the sky brightness from 14.8 mag/as2 Vega --> desired
        * level (currently 15.8 mag/as2 Vega, see UKIRT User Guide).
        */
      count *= airmass * exp(-gsAtmContOp(obs,lambda)*airmass/1.086)
                * exp((14.8-15.8)/1.086);

      /* This line is within the range of this spectrograph arm.
        * Need to get how the counts are distributed and then add this line.
        */
      iref = (long)floor(pos-SP_PSF_LEN/2+0.5);
      if (iref<0) iref=0;
      if (iref>spectro->npix[i_arm]-SP_PSF_LEN) iref=spectro->npix[i_arm]-SP_PSF_LEN;
      gsSpectroDist(spectro,i_arm,lambda,pos-iref,0,SP_PSF_LEN,FR);
      for(j=0;j<SP_PSF_LEN;j++)
        Noise[iref+j] += count*FR[j]*sample_factor;
    }
  }
}

void gsAddSkyLines(SPECTRO_ATTRIB* spectro, OBS_ATTRIB* obs, int i_arm, double* Noise) {
  printf("  --> Computing Sky Lines Contribution ...\n");
  /* Sky line contributions -- uses VLT/UVES sky model. */
  switch((obs->skytype>>16) & 0xf) {

    /* No sky lines at all, for testing only */
    case 0x0:
      break;

    case 0x1:
      gsAddSkyLines_UVES(spectro, obs, i_arm, Noise);
      gsAddSkyLines_NIR(spectro, obs, i_arm, Noise);
      break;

    default:
      gsError("Error: illegal sky line model: %1lx.", (obs->skytype>>16) & 0xf);
      break;
  }
}

double gsGetCount(SPECTRO_ATTRIB* spectro, OBS_ATTRIB* obs, int i_arm, double lambda, double continuum, double rad) {  
  double count;
  count = continuum * spectro->dl[i_arm] 
          * gsAeff(spectro,obs,i_arm)
          * gsThroughput(spectro,i_arm,lambda)
          * obs->t_exp * M_PI * rad * rad
          * gsFracTrace(spectro,i_arm,lambda,1);
  return count;
}

void gsAddLunarContinuum(SPECTRO_ATTRIB* spectro, OBS_ATTRIB* obs, int i_arm, double* Noise) {
  long ipix;
  double lambda;
  double rad, sample_factor;
  double lunar_cont, lunarphase;
  double count;
  double scale_RS, scale_MS;
  double kV;
  double Istar, f1, f2, alpha, Bmoon;

  rad = gsGetFiberRadius(spectro, obs, i_arm);
  sample_factor = gsGetSampleFactor(spectro, i_arm);

  printf("  --> Computing Lunar Continuum Contribution ...\n");
  for(ipix=0;ipix<spectro->npix[i_arm];ipix++) {
    lambda = spectro->lmin[i_arm] + (ipix+0.5)*spectro->dl[i_arm];   
    gsPrintProgress(spectro->npix[i_arm], ipix);

    /* Moonlight -- if the Moon is above the horizon */
    if (obs->lunarZA<90) {

      lunar_cont = 0.;
      lunarphase = obs->lunarphase - floor(obs->lunarphase);

      switch((obs->skytype>>4) & 0xf) {

        /* Krisciunas & Schaefer 1991 model */
        case 0x0:

          /* Color model from CFHT Redeye manual.
           *
           * Solar spectrum rescaling; right now just a 5800 K blackbody,
           * in units of photons/s/m^2/ln lambda, and then multiplied by
           * lambda^-4 (Rayleigh) or lambda^-1.3 (Mie). All of this is
           * normalized to 550 nm.
           *
           * Note: 5800 K --> kT=hc/lambda @ 2480 nm.
           */
          scale_RS = (exp(2480./550.)-1.)/(exp(2480./lambda)-1.) * pow(lambda/550., -7.0);
          scale_MS = (exp(2480./550.)-1.)/(exp(2480./lambda)-1.) * pow(lambda/550., -4.3);

          /* Model for the V band moonlight brightness */
          kV = gsAtmContOp(obs,550.);
          alpha = 360*fabs(lunarphase-0.5);
          Istar = pow(10., -0.4*(3.84 + 0.026*alpha +4e-9*alpha*alpha*alpha*alpha));
          if (alpha < 7) Istar *= 1.35-0.05*alpha;
          f1 = 2.29e5*(1.06+cos(obs->lunarangle*DEGREE)*cos(obs->lunarangle*DEGREE));
          f2 = pow(10., 6.15-obs->lunarangle/40.);
          Bmoon = (f1*scale_RS+f2*scale_MS) * Istar
                  * pow(10., -0.4*kV/sqrt(1.-0.96*sin(obs->lunarZA*DEGREE)*sin(obs->lunarZA*DEGREE)))
                  * (1. - pow(10., -0.4*kV/sqrt(1.-0.96*sin(obs->zenithangle*DEGREE)*sin(obs->zenithangle*DEGREE))));

          /* Continuum conversion at V band: [at our level of accuracy, neglect Vega-AB correction]
           * 3.408e10 = brightness in nanoLambert of a 0th magnitude object per arcsec^2
           * 5.48e10 = photons per second per m^2 per ln lambda from a 0th magnitude object
           */
          lunar_cont = Bmoon / 3.408e10 * 5.48e10 / 550.;

          break;

        default:
          gsError("Error: illegal Moonlight model: %1lx.", (obs->skytype>>4) & 0xf);
          break;
      }

    }

    count = gsGetCount(spectro, obs, i_arm, lambda, lunar_cont, rad);
    Noise[ipix] += count*sample_factor;
  }
}

double gsGetSkyContinuum_UVES(OBS_ATTRIB* obs, double lambda) {
  double continuum;

  continuum = lambda<375? 0.17: lambda<483? 0.14: lambda<580? 0.09: lambda<674.5? 0.10: lambda<858? 0.08: 0.07;
  continuum /= 1.1;
  if (lambda>1040 && ((obs->skytype & 0xf) == 0x1)) continuum = 0.4669*(1000./lambda)*( 2.0*lambda/1000. - 0.5 );
  continuum *= 1e-11 * PHOTONS_PER_ERG_1NM*lambda;

  return continuum;
}

double gsGetSkyContinuum_Gunn1(OBS_ATTRIB* obs, double lambda) {
  double mag, continuum;

  mag = 21.55 + (lambda-600)*(lambda>600?5e-5:-6e-3) - 0.55*exp(-0.005*(lambda-594)*(lambda-594))-0.175*(1.+tanh(375-lambda)) - 6.14656e9/lambda/lambda/lambda/lambda;
	/* modified by Kiyoto Yabe 20150608 */
	/*  mag = (lambda>600?24.316:27.166) + (lambda>600?-5.199e-03:-1.419e-02)*lambda + (lambda>600?1.465e-06:8.541e-06)*lambda*lambda - 0.55*exp(-0.005*(lambda-594)*(lambda-594)) - 6.14656e9/lambda/lambda/lambda/lambda;*/
  continuum = 0.01089*pow(10, 0.4*(22.5-mag)) * 1e6/lambda/lambda;
  continuum *= 1e-11 * PHOTONS_PER_ERG_1NM*lambda;

  return continuum;
}

double gsGetSkyContinuum_Gunn2(OBS_ATTRIB* obs, double lambda) {
  double mag, continuum;

  //fprintf(stderr, "Modified JG\n");
  mag = (lambda>600?24.316:27.166) + (lambda>600?-5.199e-03:-1.419e-02)*lambda + (lambda>600?1.465e-06:8.541e-06)*lambda*lambda - 0.55*exp(-0.005*(lambda-594)*(lambda-594)) - 6.14656e9/lambda/lambda/lambda/lambda;
  continuum = 0.01089*pow(10, 0.4*(22.5-mag)) * 1e6/lambda/lambda;
  continuum *= 1e-11 * PHOTONS_PER_ERG_1NM*lambda;

  return continuum;
}

double gsGetSkyContinuum_BigBoss(OBS_ATTRIB* obs, double lambda) {
  double continuum;

  continuum = 0.035*sqrt(1000/lambda) + 0.045*exp(-0.005*(lambda-594)*(lambda-594));
  continuum *= pow(10, 0.4*( 0.05+6.14656e9/lambda/lambda/lambda/lambda ));
  continuum *= 1e-11 * PHOTONS_PER_ERG_1NM*lambda;

  return continuum;
}

void gsAddSkyContinuum(SPECTRO_ATTRIB* spectro, OBS_ATTRIB* obs, int i_arm, double* Noise) {
  
  long ipix;
  double lambda;
  double airmass, rad, sample_factor;
  double trans;
  double continuum, count;
  
  airmass = gsGetAirmass(obs);
  rad = gsGetFiberRadius(spectro, obs, i_arm);
  sample_factor = gsGetSampleFactor(spectro, i_arm);

  /* Sky continuum contributions. */
  printf("  --> Computing Sky Continuum Contribution ...\n");
  for(ipix=0;ipix<spectro->npix[i_arm];ipix++) {
    lambda = spectro->lmin[i_arm] + (ipix+0.5)*spectro->dl[i_arm];   
    gsPrintProgress(spectro->npix[i_arm], ipix);
    
    /* Atmospheric transmission -- used to remap the continuum model */
    trans = gsAtmTransInst(spectro, obs, i_arm, lambda);

    /* Continuum formula -- in photons/s/m^2/arcsec^2/nm.
     * Last 4 bits of skytype determine sky model.
     */
    switch(obs->skytype & 0xf) {

      /* No sky continuum at all, for testing only */
      case 0x0:
        continuum = 0.;
        break;

      /* UVES continuum, typical airmass ~ 1.1, supplemented with IR. */
      case 0x1:
      case 0x2:
        continuum = gsGetSkyContinuum_UVES(obs, lambda);
        break;

      /* Fit to Jim Gunn spectrum */
      case 0x3:
        continuum = gsGetSkyContinuum_Gunn1(obs, lambda);
        break;

      /* Fit to Big Boss proposal spectrum */
      case 0x4:
        continuum = gsGetSkyContinuum_BigBoss(obs, lambda);
        break;

	    /* Fit to Jim Gunn spectrum (modified) */
	    /*Added by Kiyoto Yabe and modified by Yuki Moritani 20150608*/
      case 0x5:
        continuum = gsGetSkyContinuum_Gunn2(obs, lambda);
        break;

      default:
        gsError("Error: illegal sky continuum model: %1lx.", obs->skytype & 0xf);
        break;
    }

    continuum *= airmass * trans;
    count = gsGetCount(spectro, obs, i_arm, lambda, continuum, rad);
    Noise[ipix] += count*sample_factor;
  }
}

void gsAddSkySubtrError(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double* Noise, double* sky) {
  long ipix, j;
  double sky_sysref;

  printf("  --> Computing Sky Systematic Error Contribution ...\n");
  /* Add the systematic sky subtraction error */
  for(ipix=0;ipix<spectro->npix[i_arm];ipix++) {
    sky_sysref = sky[ipix];
    for(j=ipix-1;j<=ipix+1;j++)
      if (j>=0 && j<spectro->npix[i_arm])
        if (sky[j]>sky_sysref)
          sky_sysref=sky[j];
    Noise[ipix] += spectro->sysfrac*spectro->sysfrac*sky_sysref*sky_sysref;
  }
}

void gsAddStrayLight(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double* Noise, double* sky) {
  long ipix;
  double sample_factor;
  double sky_sysref;
  
  sample_factor = gsGetSampleFactor(spectro, i_arm);

  /* Add the diffuse stray light background */
  sky_sysref = 0.;
  for(ipix=0;ipix<spectro->npix[i_arm];ipix++)
    sky_sysref += sky[ipix];
  sky_sysref *= spectro->width[i_arm]*spectro->pix[i_arm]/spectro->sep[i_arm]/(double)spectro->npix[i_arm];
  for(ipix=0;ipix<spectro->npix[i_arm];ipix++)
     Noise[ipix] += spectro->diffuse_stray*sky_sysref*sample_factor;
}

void gsAddDarkNoise(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double* Noise) {
  long ipix;
  double sample_factor;
  double var;
  
  sample_factor = gsGetSampleFactor(spectro, i_arm);

  var = (spectro->dark[i_arm]*obs->t_exp*sample_factor) * spectro->width[i_arm];
  for(ipix=0;ipix<spectro->npix[i_arm];ipix++)
     Noise[ipix] += var;
}

void gsAddReadoutNoise(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double* Noise) {
  long ipix;
  double var;

  var = (spectro->read[i_arm]*spectro->read[i_arm]) * spectro->width[i_arm];
  for(ipix=0;ipix<spectro->npix[i_arm];ipix++)
     Noise[ipix] += var * obs->n_exp;
}

/* Routine to construct the noise in a given spectrograph arm.
 * Returns noise variance in counts^2 per l-pixel.
 */
void gsGetNoise(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double *Noise, double *SkyMod) {

  long ipix;
  double sample_factor;
  double *sky;

  sample_factor = gsGetSampleFactor(spectro, i_arm);

  printf(" //Arm%d//\n",i_arm);

  for(ipix=0;ipix<spectro->npix[i_arm];ipix++) SkyMod[ipix] = 0.;
  for(ipix=0;ipix<spectro->npix[i_arm];ipix++) Noise[ipix] = 0.;

  /* Sky lines */
  gsAddSkyLines(spectro, obs, i_arm, Noise);

  /* Sky and lunar continuum */
  gsAddSkyContinuum(spectro, obs, i_arm, Noise);
  gsAddLunarContinuum(spectro, obs, i_arm, Noise);

  /* Compute the sky, and add systematic error contribution */
  /* Allocation of sky vector */
  sky = (double*)malloc((size_t)(MAXPIX*sizeof(double)));
  for(ipix=0;ipix<spectro->npix[i_arm];ipix++) {
    sky[ipix]=Noise[ipix]/sample_factor;
    SkyMod[ipix]=Noise[ipix]/sample_factor;
  }

  /* Add the systematic sky subtraction error */
  gsAddSkySubtrError(spectro, obs, i_arm, Noise, sky);

  /* Add the diffuse stray light background */
  gsAddStrayLight(spectro, obs, i_arm, Noise, sky);

  /* Dark current & read noise contributions */
  gsAddDarkNoise(spectro, obs, i_arm, Noise);
  gsAddReadoutNoise(spectro, obs, i_arm, Noise);

  free((char*)sky);

  return;
}

/* Computes the signal for a feature with specified flux and velocity
 * dispersion. Output is Signal[0..Npix-1]. (Npix from spectro)
 *
 * lambda = wavelength of feature in nm.
 * F = flux in erg/cm2/s
 * sigma_v = velocity dispersion in km/s
 * r_eff = source effective radius in arcsec
 * decent (decenter of fiber from target in arcsec)
 * fieldang (field angle in degrees; 0 on axis)
 */
void gsGetSignal(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double lambda,
  double F, double sigma_v, double *Signal) {

  /* Do computations only over a finite interval, in this case 32 pixels around the feature. */
#define NP_WIN 32
  long ipix,iref;
  double pos, counts;
  double FR[NP_WIN];
  double trans, den;
  double x;

  for(ipix=0;ipix<spectro->npix[i_arm];ipix++) Signal[ipix] = 0.;

  /* Find feature location; exit if no signal */
  pos = (lambda-spectro->lmin[i_arm])/spectro->dl[i_arm];
  iref = (long)floor(pos-NP_WIN/2.0);
  if (iref<-NP_WIN || iref>=spectro->npix[i_arm]) return;
  if (iref<0) iref=0;
  if (iref>spectro->npix[i_arm]-NP_WIN) iref=spectro->npix[i_arm]-NP_WIN;

  /* Atmospheric transmission */
  trans = den = 0.;
  for(x=-4;x<4.01;x+=.2) {
    trans += gsAtmTrans(obs,lambda*(1 + x*sigma_v/299792.458)) * exp(-0.5*x*x);
    den += exp(-0.5*x*x);
  }
  trans /= den;

  /* Determine how many counts we get from the object */
  counts = F * trans
           * pow(10., -0.4*gsGalactic_Alambda__EBV(lambda)*obs->EBV)
           * gsGeometricThroughput(spectro,obs,lambda)
           * gsFracTrace(spectro,i_arm,lambda,0)
           * PHOTONS_PER_ERG_1NM * lambda * obs->t_exp 
           * gsAeff(spectro,obs,i_arm)
           * gsThroughput(spectro,i_arm,lambda) * 1e4;

  /* Get distribution of light over pixels */
  gsSpectroDist(spectro,i_arm,lambda,pos-iref,sigma_v/299792.458*lambda/spectro->dl[i_arm],NP_WIN,FR);
  for(ipix=0;ipix<NP_WIN;ipix++) Signal[ipix+iref] = FR[ipix]*counts;
  return;
}

/* Generates the signal/noise ratio for a spectral line given the noise vector Noise.
 * [This could be computed from the data given but recalculation would be expensive.]
 *
 * Line parameters are lambda (nm), F (erg/cm2/s), sigma_v (km/s), r_eff (arcsec),
 * decent (arcsec), fieldang (degrees).
 *
 * The "snrTypes" types are:
 *  0 = 1D optimal
 *  1 = uniform matched filter
 */
double gsGetSNR(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double lambda,
  double F, double sigma_v, double *Noise, int snrType) {

  long ipix,Npix;
  double SNR = 0;
  double numer, denom;
  double *Signal;

  /* Allocate and get the signal vector */
  Npix = spectro->npix[i_arm];
  Signal = (double*)malloc((size_t)(Npix*sizeof(double)));
  gsGetSignal(spectro,obs,i_arm,lambda,F,sigma_v,Signal);

  /* 1D optimal */
  if (snrType == 0) {
    SNR = 0;
    for(ipix=0;ipix<Npix;ipix++) SNR += Signal[ipix]*Signal[ipix]/Noise[ipix];
    SNR = sqrt(SNR);
  }

  /* uniform matched filter */
  if (snrType == 1) {
    numer = denom = 0;
    for(ipix=0;ipix<Npix;ipix++) {
      numer += Signal[ipix]*Signal[ipix];
      denom += Signal[ipix]*Signal[ipix]*Noise[ipix];
    }
    SNR = numer>=0.001? numer/sqrt(denom): 0;
  }

  free((char*)Signal);
  return(SNR);
}

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

/* Get atmospheric transmission for the continuum by taking instrumental effects into account.
 */
double gsAtmTransInst(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double lambda) {
  int j;
  double FR[5*SP_PSF_LEN];
  double den, num;
  double trans;

  /* Atmospheric transmission */
  gsSpectroDist(spectro,i_arm,lambda,7.5,0,SP_PSF_LEN,FR);
  num = den = 0.0;
  for(j=0;j<5*SP_PSF_LEN;j++) {
    trans = gsAtmTrans(obs,lambda+(0.2*j-SP_PSF_LEN/2+0.5)*spectro->dl[i_arm]);
    num += FR[j/5]*trans;
    den += FR[j/5];
  }
  trans = num/den;
  return trans;
}

double gsMagToFlux(double mag) {
  double flux;
  flux = 3.631e-20 * pow(10., -0.4*mag); /* in erg/cm2/s */
  return flux;
}

double gsEBVToAttn(OBS_ATTRIB* obs, double lambda) {
  double attn;    /* attenuation at lambda */
  attn = pow(10., -0.4*gsGalactic_Alambda__EBV(lambda)*obs->EBV);
  return attn;
}

double gsPerHertzToPerPixel(SPECTRO_ATTRIB* spectro, int i_arm, double lambda) {
  return 2.99792458e17*spectro->dl[i_arm]/(lambda*lambda);
}

/* Compute overall transmission function including reddening, atmosphere and instrument
 * Do conversion between erg/cm2/s to counts per pixel 
 * r_eff is galaxy effective radius [arcsec], set 0 for stars
 * decent is fiber astrometric offset [arcsec]
*/
double gsConversionFunction(SPECTRO_ATTRIB *spectro, OBS_ATTRIB *obs, int i_arm, double lambda) {
  double conv;
  
  /* Atmospheric transmission */
  conv = gsAtmTransInst(spectro, obs, i_arm, lambda);
  /* Galactic reddening */
  conv *= gsEBVToAttn(obs, lambda);
  /* Instrument geometry */
  conv *= gsGeometricThroughput(spectro,obs,lambda);
  /* Line spread function */
  conv *= gsFracTrace(spectro,i_arm,lambda,0);
  /* Effective area, including throughput, and time */
  conv *= PHOTONS_PER_ERG_1NM * lambda * obs->t_exp 
          * gsAeff(spectro,obs,i_arm)
          * gsThroughput(spectro,i_arm,lambda) * 1e4;
  /* Convert from per Hz --> l-per pixel */
  conv *= gsPerHertzToPerPixel(spectro, i_arm, lambda);
  return conv;
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

/* --- I/O FUNCTIONS --- */

int spectro_arm(const SPECTRO_ATTRIB *spectro, int ia)
{
   if (!spectro->MR) {
      return ia;
   } else {
      return (ia == 1) ? 3 : ia;
   }
}

void gsPrintCompilerFlags() {
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
}

/* Open a config file and quit is does not exist.
 * If file name is -, open stdin
 */
FILE* gsOpenConfigFile(const char* filename) {
  if (strcmp("-", filename) == 0) {
    return stdin;
  } else {
    FILE* fp;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        gsError("Error: Can't read file: %s.", filename);
    }
    return fp;
  }
}

void gsCloseConfigFile(FILE* fp) {
  fclose(fp);
}

FILE* gsOpenOutputFile(const char* filename) {
  if (strcmp("-", filename) == 0) {
    return stdout;
  } else {
    FILE* fp;
    fp = fopen(filename, "w");
    if (fp == NULL) {
        gsError("Error: Can't open file for writing: %s.", filename);
    }
    return fp;
  }
}

void gsCloseOutputFile(FILE* fp) {
  fclose(fp);
}

void gsReadObservationKeyword(char InfoLine[512], const char* keyword, const char* format, void* var)
{
  int args;
  int len = strlen(keyword);
  if (memcmp(InfoLine, keyword, (size_t)len)==0) {
    args = sscanf(InfoLine+len+1, format, var);
    if (args!=1) {
      gsError("Error: gsReadObservationConfig: Failed to read %s keyword: %d/1 arguments assigned.", keyword, args);
    }
  }
}

/* Read an observation configuration file
 */
void gsReadObservationConfig(FILE *fp, OBS_ATTRIB *obs) {
  char InfoLine[512];

  obs->seeing_fwhm_800 = 0.8;
  obs->elevation = 4000.0;
  obs->zenithangle = 45.0;
  obs->lunarphase = 0.3;
  obs->lunarangle = 120.0;
  obs->lunarZA = 70.0;
  obs->EBV = 0.4;
  obs->fieldangle = 0.675;
  obs->decent = 0;
  obs->t_exp = 900.0;
  obs->n_exp = 4;
  obs->r_eff = 0.0;
  obs->skytype = 0x00011005;
  obs->flags = 0;

  do {
    if (fgets(InfoLine, 511, fp)==NULL) {
      InfoLine[0] = '@';
    }
    else {
      gsReadObservationKeyword(InfoLine, "SEEING", "%lg", &(obs->seeing_fwhm_800));
      gsReadObservationKeyword(InfoLine, "ELEV", "%lg", &(obs->elevation));
      gsReadObservationKeyword(InfoLine, "ZA", "%lg", &(obs->zenithangle));
      gsReadObservationKeyword(InfoLine, "LUNARPHASE", "%lg", &(obs->lunarphase));
      gsReadObservationKeyword(InfoLine, "LUNARANGLE", "%lg", &(obs->lunarangle));
      gsReadObservationKeyword(InfoLine, "LUNARZA", "%lg", &(obs->lunarZA));
      gsReadObservationKeyword(InfoLine, "EBV", "%lg", &(obs->EBV));
      gsReadObservationKeyword(InfoLine, "FIELDANG", "%lg", &(obs->fieldangle));
      gsReadObservationKeyword(InfoLine, "DECENT", "%lg", &(obs->decent));
      gsReadObservationKeyword(InfoLine, "TEXP", "%lg", &(obs->t_exp));
      gsReadObservationKeyword(InfoLine, "NEXP", "%d", &(obs->n_exp));
      gsReadObservationKeyword(InfoLine, "REFF", "%lg", &(obs->r_eff));
    }
  } while (InfoLine[0]!='@');
}

/* Write a single observation config file keyword with optional prefix
 */
void gsWriteObservationKeyword(FILE* fp, const char* prefix, const char* keyword, const char* format, ...) {
  if (prefix != NULL) 
    fprintf(fp, "%s", prefix);
  fprintf(fp, "%s ", keyword);

  va_list argptr;
  va_start(argptr, format);
  vfprintf(fp, format, argptr);
  va_end(argptr);

  fprintf(fp, "\n");
}

/* Write observation configuration file, with optional line prefix
 */
void gsWriteObservationConfig(FILE *fp, OBS_ATTRIB* obs, const char* prefix) {
  gsWriteObservationKeyword(fp, prefix, "SEEING", "%lg", obs->seeing_fwhm_800);
  gsWriteObservationKeyword(fp, prefix, "ELEV", "%lg", obs->elevation);
  gsWriteObservationKeyword(fp, prefix, "ZA", "%lg", obs->zenithangle);
  gsWriteObservationKeyword(fp, prefix, "LUNARPHASE", "%lg", obs->lunarphase);
  gsWriteObservationKeyword(fp, prefix, "LUNARANGLE", "%lg", obs->lunarangle);
  gsWriteObservationKeyword(fp, prefix, "LUNARZA", "%lg", obs->lunarZA);
  gsWriteObservationKeyword(fp, prefix, "EBV", "%lg", obs->EBV);
  gsWriteObservationKeyword(fp, prefix, "FIELDANG", "%lg", obs->fieldangle);
  gsWriteObservationKeyword(fp, prefix, "DECENT", "%lg", obs->decent);
  gsWriteObservationKeyword(fp, prefix, "TEXP", "%lg", obs->t_exp);
  gsWriteObservationKeyword(fp, prefix, "NEXP", "%d", obs->n_exp);
  gsWriteObservationKeyword(fp, prefix, "REFF", "%lg", obs->r_eff);
}

void gsCopyObservationConfig(OBS_ATTRIB* obs1, OBS_ATTRIB* obs2) {
  memcpy(obs2, obs1, sizeof(OBS_ATTRIB));
}

/* Read a spectrograph configuration file.
 * Format:
 * # all lines starting with # are ignored.
 * # all other lines start with a keyword.
 * OPTICS <primary diam> <central obscuration> <EFL> <FOV rad/deg>
 * SPOT <rms/axis @ center> <> <> <> <rms/axis @ edge>
 * VIGNET <center> <> <> <> <> <edge> [this line is optional; if excluded 1 is assumed]
 * FIBER <entry radius/um>
 * ARMS <N arms>
 * PARAM <arm #> <lmin> <lmax> <npix> <width>
 * CAMERA <arm #> <f/ratio> <ccd thickness> <pix scale> <T> <rms spot> <fiber image diam> <dark> <read> <trace sep>
 * THRPUT <N thr>
 * <lambda_0> <Thr(lambda_0) [1st]> ... <Thr(lambda_0) [5th]>
 * ...
 * <lambda_N-1> <Thr(lambda_N-1) [1st]> ... <Thr(lambda_N-1) [5th]>
 */
void gsReadSpectrographConfig(FILE *fp, SPECTRO_ATTRIB *spectro, double degrade) {
  int i, i_arm, args;
  char InfoLine[512];
  double *ptr, lmin, lmax;
  long npix, width, count;
  double f, thick, pix, T, rms, diam, dark, read, sep, temp[5];

  /* Set parameters to defaults or illegal values */
  spectro->D_outer = -1;
  for(i=0;i<5;i++) spectro->rms_spot[i] = -1;
  spectro->fiber_ent_rad = -1;
  spectro->N_arms = 0;
  for(i=0;i<5;i++) spectro->vignette[i] = 1.;
  for(i=0;i<MAXARM;i++) spectro->Dtype[i] = 0;

  do {
    if (fgets(InfoLine, 511, fp)==NULL) {
      InfoLine[0] = '@';
    }
    else {
      if (memcmp(InfoLine, "OPTICS", (size_t)6)==0) {
        args = sscanf(InfoLine+7, "%lg %lg %lg %lg %lg %lg %lg %lg", &(spectro->D_outer), &(spectro->centobs), &(spectro->rfov),
          spectro->EFL, spectro->EFL+1, spectro->EFL+2, spectro->EFL+3, spectro->EFL+4);
        if (args!=8) {
          gsError("Error: gsReadSpectrographConfig: Failed to read OPTICS keyword: %d/4 arguments assigned.", args);
        }
      }
      if (memcmp(InfoLine, "SPOT", (size_t)4)==0) {
        ptr = spectro->rms_spot;
        args = sscanf(InfoLine+5, "%lg %lg %lg %lg %lg", ptr, ptr+1, ptr+2, ptr+3, ptr+4);
        if (args!=5) {
          gsError("Error: gsReadSpectrographConfig: Failed to read SPOT keyword: %d/5 arguments assigned.", args);
        }
      }      
      if (memcmp(InfoLine, "FIBER", (size_t)5)==0) {
        args = sscanf(InfoLine+6, "%lg", &(spectro->fiber_ent_rad));
        if (args!=1) {
          gsError("Error: gsReadSpectrographConfig: Failed to read FIBER keyword: %d/1 arguments assigned.", args);
        }
      }      
      if (memcmp(InfoLine, "ARMS", (size_t)4)==0) {
        args = sscanf(InfoLine+5, "%d", &(spectro->N_arms));
        if (args!=1) {
          gsError("Error: gsReadSpectrographConfig: Failed to read ARMS keyword: %d/1 arguments assigned.", args);
        }
        for(i=0;i<spectro->N_arms;i++) {
          spectro->lmin[i] = spectro->fratio[i] = -1.;
          spectro->nline[i] = 1e12;
            /* Effective number of lines might as well be infinite; this is changed later by NLINES keyword */
        }
      }
      if (memcmp(InfoLine, "MEDIUM_RESOLUTION", (size_t)17)==0) {
        args = sscanf(InfoLine+17, "%d", &spectro->MR);
        if (args!=1) {
            gsError("Error: gsReadSpectrographConfig: Failed to read MEDIUM_RESOLUTION keyword.");
        } else if(spectro->MR != 0 && spectro->MR != 1) {
            gsError("Error: gsReadSpectrographConfig: invalid value of MEDIUM_RESOLUTION %d (use 0/1).", spectro->MR);
    	  }
      }
      if (memcmp(InfoLine, "PARAM", (size_t)5)==0) {
        args = sscanf(InfoLine+6, "%d %lg %lg %ld %ld", &i, &lmin, &lmax, &npix, &width);
        if (args!=5) {
          gsError("Error: gsReadSpectrographConfig: Failed to read PARAM keyword: %d/5 arguments assigned.", args);
        }
        if (i<0 || i>=spectro->N_arms) {
          gsError("Error: illegal PARAM line: arm #%1d does not exist.", i);
        }
        spectro->lmin[i] = lmin;
        spectro->lmax[i] = lmax;
        spectro->npix[i] = npix;
        spectro->dl[i] = (lmax-lmin)/(double)npix;
        spectro->width[i] = width;
      }
      if (memcmp(InfoLine, "CAMERA", (size_t)6)==0) {
        args = sscanf(InfoLine+7, "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg", &i, &f, &thick, &pix, &T, &rms, &diam, &dark, &read, &sep);
        if (args!=10) {
          gsError("Error: gsReadSpectrographConfig: Failed to read CAMERA keyword: %d/10 arguments assigned.", args);
        }
        if (i<0 || i>=spectro->N_arms) {
          gsError("Error: illegal CAMERA line: arm #%1d does not exist.", i);
        }
        spectro->fratio[i] = f;
        spectro->thick[i] = thick;
        spectro->pix[i] = pix;
        spectro->temperature[i] = T;
        spectro->rms_cam[i] = rms;
        spectro->diam[i] = diam;
        spectro->dark[i] = dark;
        spectro->read[i] = read;
        spectro->sep[i] = sep;
      }
      if (memcmp(InfoLine, "THRPUT", (size_t)6)==0) {
        i = 0;
        i_arm = 0;
        spectro->istart[0] = 0;
        for(count=0;count<1000000;count++) {
          if (fgets(InfoLine, 511, fp)==NULL) {
            gsError("Error: Unexpected EOF: @ THRPUT grid point %d.", i);
          }
          if (memcmp(InfoLine, "D", (size_t)1)==0) {
            i_arm++;
            spectro->istart[i_arm] = i;
            if (i_arm == spectro->N_arms) {
              spectro->N_thr = i;
              break;
            }
            continue;
          }
          args = sscanf(InfoLine, "%lg %lg %lg %lg %lg %lg", spectro->l+i, temp, temp+1, temp+2, temp+3, temp+4);
          if (args!=6) {
            gsError("Error: gsReadSpectrographConfig: illegal throughput table line: %d/6 arguments assigned.", args);
          }
          spectro->T[i] = temp[0]*temp[1]*temp[2]*temp[3]*temp[4]*degrade;
          i++;
        }
        if (count==1000000) {
          gsError("Error: No end of throughput table found after 1000000 lines.");
        }
      }
      if (memcmp(InfoLine, "VIGNET", (size_t)6)==0) {
        ptr = spectro->vignette;
        args = sscanf(InfoLine+7, "%lg %lg %lg %lg %lg", ptr, ptr+1, ptr+2, ptr+3, ptr+4);
        if (args!=5) {
          gsError("Error: gsReadSpectrographConfig: Failed to read VIGNET keyword: %d/5 arguments assigned.", args);
        }
      }
      if (memcmp(InfoLine, "HGCDTE", (size_t)6)==0) {
        sscanf(InfoLine+7, "%d", &i);
        if (i<0 || i>=spectro->N_arms) {
          gsError("Error: HGCDTE %d: illegal arm index.", i);
        }
        spectro->Dtype[i] = 1;
      }
      if (memcmp(InfoLine, "NLINES", (size_t)6)==0) {
        sscanf(InfoLine+7, "%d %lg", &i, temp);
        if (i<0 || i>=spectro->N_arms) {
          gsError("Error: NLINES %d: illegal arm index.", i);
        }
        if (*temp<10) {
          gsError("Error: NLINES %d: nlines=%lg is illegal.", i, *temp);
        }
        spectro->nline[i] = *temp;
      }
    }
  } while (InfoLine[0]!='@');

  /* Tests */
  if (spectro->D_outer<=0) {
    gsError("Error: illegal outer diameter or no OPTICS line.");
  }
  for(i=0; i<5; i++) if (spectro->rms_spot[i]<0) {
    gsError("Error: illegal rms spot or no SPOT line: spot[%1d]=%12.5lE.", i, spectro->rms_spot[i]);
  }
  if (spectro->fiber_ent_rad<=0) {
    gsError("Error: illegal fiber radius or no FIBER line.");
  }
  if (spectro->N_arms<=0) {
    gsError("Error: %d arms illegal or no ARMS line.", spectro->N_arms);
  }
  for(i=0;i<spectro->N_arms;i++) if (spectro->lmin[i]<0) {
    gsError("Error: illegal lambda min = %lg or no PARAM %1d line.", spectro->lmin[i], i);
  }

#ifdef NATURAL_NLINES
  /* Forces N_lines to be the actual number across the pupil; overrides NLINES keyword. */
  for(i=0;i<spectro->N_arms;i++) {
    spectro->nline[i] = spectro->pix[i]/spectro->fratio[i]/spectro->dl[i]*1000.;
      /* factor of 1000 needed since pixel scale is in um and wavelength in nm */
  }
#endif

  /* Initially set no systematics; override later with user input */
  spectro->sysfrac=0.;
  return;
}

void gsAllocArmVectors(SPECTRO_ATTRIB* spectro, double*** spec) {
  int ia;

  *spec = (double**)malloc((size_t)(spectro->N_arms*sizeof(double*)));
  for(ia=0; ia<spectro->N_arms; ia++) 
    (*spec)[ia] = (double*)calloc((size_t)spectro->npix[ia], sizeof(double));
}

void gsFreeArmVectors(SPECTRO_ATTRIB* spectro, double** spec) {
  int ia;

  for(ia=0; ia<spectro->N_arms; ia++) free((char*)(spec[ia]));
  free((char*)spec);
}

void gsReadMagfile(MAGFILE* magfile, char* filename) {
  FILE *in_pipe, *fp;
  int ia;
  char command[256], buf[256], dmy[32];

  /* This is an ugly linux-only hack to count number of lines */
  sprintf(command, "wc %s", filename);
  if ((in_pipe = popen(command,"r")) == NULL){
    exit(1);
  }
  if(fgets(buf,128, in_pipe)!=NULL)
    sscanf(buf, "%d %s", &magfile->num_inmag, dmy);

  magfile->lambda_inmag = (double*) malloc(sizeof(double) * magfile->num_inmag);
  magfile->mag_inmag = (double*) malloc(sizeof(double) * magfile->num_inmag);
  ia=0;
  fp=fopen(filename,"r");
  while(fgets(buf, 256, fp) != NULL){
    if(strncmp(buf,"#",1) == 0) continue;
    sscanf(buf, "%lf %lf", magfile->lambda_inmag+ia, magfile->mag_inmag+ia);
    ia++;
  }
  magfile->num_inmag = ia;
  fclose(fp);  
}

void gsCopyMagfile(MAGFILE* magfile, MAGFILE* magfile2) {
  int ia;

  /* Modified by K.Yabe 20160703 */
  /* Rewritten as a function by L.Dobos on 20191219 */
  magfile2->lambda_inmag = (double*) malloc(sizeof(double) * (magfile->num_inmag+2));
  magfile2->mag_inmag = (double*) malloc(sizeof(double) * (magfile->num_inmag+2));
  if (magfile->lambda_inmag[0]>300.0){
    magfile2->lambda_inmag[0] = 300.0;
    magfile2->mag_inmag[0] = 99.9;
  }
  else {
    magfile2->lambda_inmag[0] = magfile->lambda_inmag[0]-1.0;
    magfile2->mag_inmag[0] = magfile->mag_inmag[0];
  }
  ia=0;
  while(ia<magfile->num_inmag){
    magfile2->lambda_inmag[ia+1] = magfile->lambda_inmag[ia];
    magfile2->mag_inmag[ia+1] = magfile->mag_inmag[ia];
    if (magfile->mag_inmag[ia]<=0.0){
      magfile2->mag_inmag[ia+1] = 99.9;
    }
    else {
      magfile2->mag_inmag[ia+1] = magfile->mag_inmag[ia];
    }
    ia++;
  }
  if (magfile->lambda_inmag[ia-1]<1300.0){
    magfile2->lambda_inmag[ia+1] = 1300.0;
    magfile2->mag_inmag[ia+1] = 99.9;
  }
  else{
    magfile2->lambda_inmag[ia+1] = magfile->lambda_inmag[ia-1]+1.0;
    magfile2->mag_inmag[ia+1] = magfile->mag_inmag[ia-1];
  }
  magfile2->num_inmag = ia + 2;
  /* printf("%d %d\n", magfile->num_inmag, magfile2->num_inmag2); */
}

void gsFreeMagfile(MAGFILE* magfile) {
  free(magfile->lambda_inmag);
  free(magfile->mag_inmag);
}

double gsInterpolateMagfile(MAGFILE* magfile, double lambda) {
  int p1, p2, k, kk;
  double mag;
  
  /* Added by Y.Moritani for input mag. file: 20150422 */
  /* interpolate magnitude for a given lambda */
  /* mag is used as a fllag: -99.9 ... use input file*/
  /* Updated by to be a function and use a struct as input L.Dobos on 20191219 */
  kk=0;
  if(lambda < magfile->lambda_inmag[0]){
    p1=0; p2=1;
  } else if(lambda > magfile->lambda_inmag[magfile->num_inmag-1]){
    p1=magfile->num_inmag-2; p2=magfile->num_inmag-1;
  } else{
    for(k=kk; k<magfile->num_inmag-1; k++){
      if(lambda > magfile->lambda_inmag[k] && lambda <= magfile->lambda_inmag[k+1]){
        p1=k ; p2 = k+1;
        kk=k;
        //fprintf(stderr, "%lf %lf\n",magfile2->lambda_inmag[k], magfile2->lambda_inmag[k+1]);
      }
    } 
  }
  mag = ((lambda-magfile->lambda_inmag[p1])*magfile->mag_inmag[p2] + (magfile->lambda_inmag[p2]-lambda)*magfile->mag_inmag[p1])/(magfile->lambda_inmag[p2]-magfile->lambda_inmag[p1]);
  return mag;
}

void gsReadObsConfig_Legacy(OBS_ATTRIB* obs, SPECTRO_ATTRIB* spectro) {
  /* Made into a function by L.Dobos on 20191219 */
  /* Input observing conditions from stdin */

    /* Get observational conditions */
  /* Modified by Y.Moritnani -- 2016.02.16 */

  //printf("Enter observational conditions [hexadecimal code; suggested=11005]: ");
  if(scanf("%lx", &(obs->skytype))==EOF)
    obs->skytype = 0x10005;
  /*
  #if 0
    obs.skytype = 0x10005;
  #endif
  */
  /* Input observing conditions */

  //printf("Enter seeing [arcsec FWHM @ lambda=800nm]: ");
  if(scanf("%lg", &(obs->seeing_fwhm_800))==EOF)
    obs->seeing_fwhm_800 = 0.8; 

  //printf("Enter zenith angle [degrees]: ");
  if(scanf("%lg", &(obs->zenithangle))==EOF)
    obs->zenithangle=45.;

  //printf("Enter Galactic dust column [magnitudes E(B-V)]: ");
  if(scanf("%lg", &(obs->EBV))==EOF)
    obs->EBV = 0.00;

  //printf("Enter field angle [degrees]: ");
  if(scanf("%lg", &(obs->fieldangle))==EOF)
    obs->fieldangle = 0.675;

  //printf("Enter fiber astrometric offset [arcsec]: ");
  if(scanf("%lg", &(obs->decent))==EOF)
    obs->decent = 0.00;

  /* Moon conditions -- set Moon below horizon unless otherwise indicated */
  obs->lunarZA = 135.;
  obs->lunarangle = 90.;
  obs->lunarphase = 0.25;
#ifdef MOONLIGHT_
  //printf("Enter Moonlight conditions:\n");
  //printf("  Angle (Moon-Zenith) [deg]: ");
  if(scanf("%lg", &(obs->lunarZA))==EOF)
    obs->lunarZA=30.;
  //printf("  Angle (Moon-target) [deg]: ");
  if(scanf("%lg", &(obs->lunarangle))==EOF)
    obs->lunarangle=60.;
  //printf("  Lunar phase [0=New, 0.25=quarter, 0.5=full]: ");
  if(scanf("%lg", &(obs->lunarphase))==EOF)
    obs->lunarphase=0.;
#endif
    
  /* Exposure time and systematics */
  //printf("Enter time per exposure [s]: ");
  if(scanf("%lg", &(obs->t_exp))==EOF)
    obs->t_exp=450.;
  //printf("Enter number of exposures: ");
  if(scanf("%d", &(obs->n_exp))==EOF)
    obs->n_exp=8;
  //printf("Enter systematic sky subtraction floor [rms per 1D pixel]: ");
  if(scanf("%lg", &(spectro->sysfrac))==EOF)
    spectro->sysfrac=0.01;
  spectro->sysfrac *= sqrt((double)obs->n_exp); /* Prevent this from averaging down */
  //printf("Enter diffuse stray light [fraction of total]: ");
  if(scanf("%lg", &(spectro->diffuse_stray))==EOF)
    spectro->diffuse_stray=0.02;

  obs->flags = 0x0;
}

/* Command-line argument parsing */
char* gsGetArgPositional(int argc, char* argv[], int i) {
  if (i < 1 || argc <= i || (argv[i][0] == '-' && argv[i][1] == '-')) {
    gsError("Positional argument %d must be specified.", i);
  }
  return argv[i];
}

int gsGetArgNamedBoolean(int argc, char* argv[], const char* name) {
  int i;
  for (i = 1; i < argc; i++) {
    if (strcmp(name, argv[i]) == 0) {
      return 1;
    }
  }
  return 0;
}

void gsError(const char* message, ...) {
    va_list argptr;
    va_start(argptr, message);
    vfprintf(stderr, message, argptr);
    va_end(argptr);
    fprintf(stderr, "\n");
    exit(1);
}

void gsPrintProgress(long total, long i) {
#ifndef NOPROGRESS
  printf("      --> %.0f percent done ...\r", (double)(i + 1) / total * 100);
#endif
}