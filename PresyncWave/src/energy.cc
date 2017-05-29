#include "cctk.h" 
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void presync_energy(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
      
  const int imin0=cctk_nghostzones[0];
  const int imin1=cctk_nghostzones[1];
  const int imin2=cctk_nghostzones[2];
  const int imax0=cctk_lsh[0] - cctk_nghostzones[0];
  const int imax1=cctk_lsh[1] - cctk_nghostzones[1];
  const int imax2=cctk_lsh[2] - cctk_nghostzones[2];
  const int zero = CCTK_GFINDEX3D(cctkGH,0,0,0);
  const int di = (cctk_lsh[0]==1) ? 0:CCTK_GFINDEX3D(cctkGH,1,0,0) - zero;
  const int dj = (cctk_lsh[1]==1) ? 0:CCTK_GFINDEX3D(cctkGH,0,1,0) - zero;
  const int dk = (cctk_lsh[2]==1) ? 0:CCTK_GFINDEX3D(cctkGH,0,0,1) - zero;
  #pragma omp parallel
  CCTK_LOOP3(calc_energy,i,j,k,
    imin0,imin1,imin2,imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    int cc = CCTK_GFINDEX3D(cctkGH,i,j,k);
    double psix = (psi[cc+di]-psi[cc-di])/(2.0*CCTK_DELTA_SPACE(0));
    double psiy = (psi[cc+dj]-psi[cc-dj])/(2.0*CCTK_DELTA_SPACE(1));
    double psiz = (psi[cc+dk]-psi[cc-dk])/(2.0*CCTK_DELTA_SPACE(2));
    energy[cc] = phi[cc]*phi[cc] + psix*psix +
      psiy*psiy + psiz*psiz;
  }
  CCTK_ENDLOOP3(calc_energy);
  
}
