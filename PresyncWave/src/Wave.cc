#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <iostream>

extern "C" void Boundary_ApplyPhysicalBCs(CCTK_ARGUMENTS);

extern "C" void Carpet_ApplyPhysicalBCs(const cGH *cctkGH);

#define sq(X) (X)*(X)

extern "C"
void ready_step(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {
        int cc = CCTK_GFINDEX3D(cctkGH,i,j,k);
        //psi[cc] = psi_p[cc];
        //phi[cc] = phi_p[cc];
      }
    }
  }
}

extern "C"
void presync_wave_init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  const int imin0=0;//cctk_nghostzones[0];
  const int imin1=0;//cctk_nghostzones[1];
  const int imin2=0;//cctk_nghostzones[2];
  const int imax0=cctk_lsh[0];// - cctk_nghostzones[0];
  const int imax1=cctk_lsh[1];// - cctk_nghostzones[1];
  const int imax2=cctk_lsh[2];// - cctk_nghostzones[2];
  const int zero = CCTK_GFINDEX3D(cctkGH,0,0,0);
  const int di = CCTK_GFINDEX3D(cctkGH,1,0,0) - zero;
  const int dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - zero;
  const int dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - zero;
  CCTK_REAL x0 = x[CCTK_GFINDEX3D(cctkGH,cctk_lsh[0]/2,cctk_lsh[1]/2,cctk_lsh[2]/2)];
  CCTK_REAL y0 = x[CCTK_GFINDEX3D(cctkGH,cctk_lsh[0]/2,cctk_lsh[1]/2,cctk_lsh[2]/2)];
  #pragma omp parallel
  CCTK_LOOP3(calc_presync_wave_init,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    int cc = CCTK_GFINDEX3D(cctkGH,i,j,k);
    phi[cc] = exp(-sq(x[cc]-x0)-sq(y[cc]-y0));
    psi[cc] = 0;
  }
  CCTK_ENDLOOP3(calc_presync_wave_init);
}

extern "C"
void presync_wave_evolve(CCTK_ARGUMENTS)
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
  const int di = CCTK_GFINDEX3D(cctkGH,1,0,0) - zero;
  const int dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - zero;
  const int dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - zero;
  #pragma omp parallel
  CCTK_LOOP3(calc_presync_wave_evol,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    int cc = CCTK_GFINDEX3D(cctkGH,i,j,k);
    psi_rhs[cc] = phi[cc];
    phi_rhs[cc] = dxx_psi[cc]+dyy_psi[cc]+dzz_psi[cc];
  }
  CCTK_ENDLOOP3(calc_presync_wave_evol);
}

extern "C"
void presync_derivatives(CCTK_ARGUMENTS)
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
  const int di = CCTK_GFINDEX3D(cctkGH,1,0,0) - zero;
  const int dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - zero;
  const int dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - zero;
  assert(!std::isnan(psi[zero]));
  #pragma omp parallel
  CCTK_LOOP3(calc_presync_derivs,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    int cc = CCTK_GFINDEX3D(cctkGH,i,j,k);
    dxx_psi[cc] = (psi[cc+di]+psi[cc-di]-2.0*psi[cc])/(2.0*CCTK_DELTA_SPACE(0));
    dyy_psi[cc] = (psi[cc+dj]+psi[cc-dj]-2.0*psi[cc])/(2.0*CCTK_DELTA_SPACE(1));
    dzz_psi[cc] = (psi[cc+dk]+psi[cc-dk]-2.0*psi[cc])/(2.0*CCTK_DELTA_SPACE(2));
  }
  CCTK_ENDLOOP3(calc_presync_derivs);
}
