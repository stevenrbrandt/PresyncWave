#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <algorithm>
#include "PSWave.h"

CCTK_INT fun_stwave(
  const cGH *cctkGH,
  int num_vars,
  int *var_indices,
  int *faces,
  int *widths,
  int *table_handles) {
  _DECLARE_CCTK_ARGUMENTS;

  for(int vi=0;vi<num_vars;vi++) {
    assert(widths[vi] > 0);
    int var_index = var_indices[vi];
    for(int face=0;face<6;face++) {

      if(!cctk_bbox[face])
        continue; // we aren't on a physical boundary

/* These statements check for trivial dimensions with only one grid point.
   While the code is full 3D, it can handle the special case of setting one or two
   dimensions to one grid point as long as the initial conditions and par file are
   valid. */
      int faceValues[] = {0,1,2,3,4,5};
      if (( std::find(faceValues,faceValues+2,face) != faceValues+2 ) && cctk_lsh[0] == 1) 
        continue;

      if (( std::find(faceValues+2,faceValues+4,face) != faceValues+4 ) && cctk_lsh[1] == 1) 
        continue;
      
      if (( std::find(faceValues+4,faceValues+6,face) != faceValues+6 ) && cctk_lsh[2] == 1) 
        continue;

      int width = widths[vi];
      assert(width > 0);
      CCTK_REAL *var = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH,0,var_index);
     if(face == 0) {
        for(int k=0;k<cctk_lsh[2];k++) {
          for(int j=0;j<cctk_lsh[1];j++) {
            int cc0 = CCTK_GFINDEX3D(cctkGH,0,j,k);
            var[cc0] = 0;
          }
        }
      } else if(face == 1) {
        for(int k=0;k<cctk_lsh[2];k++) {
          for(int j=0;j<cctk_lsh[1];j++) {
            int cc0 = CCTK_GFINDEX3D(cctkGH,cctk_lsh[0]-1,j,k);
            var[cc0] = 0;
          }
        }
      } else if(face == 2) {
        for(int k=0;k<cctk_lsh[2];k++) {
          for(int i=0;i<cctk_lsh[0];i++) {
            int cc0 = CCTK_GFINDEX3D(cctkGH,i,0,k);
            var[cc0] = 0;
          }
        }
      } else if(face == 3) {
        for(int k=0;k<cctk_lsh[2];k++) {
          for(int i=0;i<cctk_lsh[0];i++) {
            int cc0 = CCTK_GFINDEX3D(cctkGH,i,cctk_lsh[1]-1,k);
            var[cc0] = 0;
          }
        }
      } else if(face == 4) {
        for(int j=0;j<cctk_lsh[1];j++) {
          for(int i=0;i<cctk_lsh[0];i++) {
            int cc0 = CCTK_GFINDEX3D(cctkGH,i,j,0);
            var[cc0] = 0;
          }
        }
      } else if(face == 5) {
        for(int j=0;j<cctk_lsh[1];j++) {
          for(int i=0;i<cctk_lsh[0];i++) {
            int cc0 = CCTK_GFINDEX3D(cctkGH,i,j,cctk_lsh[2]-1);
            var[cc0] = 0;
          }
        }
      } else {
        abort();
      }
    }
  }
  return 0;
}

// 1 = symmetric, -1 = antisymmetric
// 1 = straddle: point is a grid point
// 0 = straddle: point is between grid points.

template<int sym,int straddle>
CCTK_INT fun_bf2(
  const cGH *cctkGH,
  int num_vars,
  int *var_indices,
  int *faces,
  int *widths,
  int *table_handles) {
  _DECLARE_CCTK_ARGUMENTS;
  for(int vi=0;vi<num_vars;vi++) {
    assert(widths[vi] > 0);
    int var_index = var_indices[vi];
    for(int face=0;face<6;face++) {

      if(!cctk_bbox[face])
        continue; // we aren't on a physical boundary

/* These statements check for trivial dimensions with only one grid point.
   While the code is full 3D, it can handle the special case of setting one or two
   dimensions to one grid point as long as the initial conditions and par file are
   valid. */
      int faceValues[] = {0,1,2,3,4,5};
      if (( std::find(faceValues,faceValues+2,face) != faceValues+2 ) && cctk_lsh[0] == 1) 
        continue;

      if (( std::find(faceValues+2,faceValues+4,face) != faceValues+4 ) && cctk_lsh[1] == 1) 
        continue;
      
      if (( std::find(faceValues+4,faceValues+6,face) != faceValues+6 ) && cctk_lsh[2] == 1) 
        continue;

      int width = widths[vi];
      assert(width > 0);
      CCTK_REAL *var = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH,0,var_index);
      if(face == 0) {
        for(int k=0;k<cctk_lsh[2];k++) {
          for(int j=0;j<cctk_lsh[1];j++) {
            for(int i=0;i<width;i++) {
              int i0 = i;
              int i2 = 2*width-i+1-straddle;
              int cc0 = CCTK_GFINDEX3D(cctkGH,i0,j,k);
              int cc2 = CCTK_GFINDEX3D(cctkGH,i2,j,k);
              var[cc0] = sym*var[cc2];
            }
          }
        }
      } else if(face == 1) {
        for(int k=0;k<cctk_lsh[2];k++) {
          for(int j=0;j<cctk_lsh[1];j++) {
            for(int i=0;i<width;i++) {
              int i0 = cctk_lsh[0]-1-i;
              int i2 = cctk_lsh[0]-2*width+i-2+straddle;
              int cc0 = CCTK_GFINDEX3D(cctkGH,i0,j,k);
              int cc2 = CCTK_GFINDEX3D(cctkGH,i2,j,k);
              var[cc0] = sym*var[cc2];
            }
          }
        }
      } else if(face == 2) {
        for(int k=0;k<cctk_lsh[2];k++) {
          for(int j=0;j<width;j++) {
            for(int i=0;i<cctk_lsh[0];i++) {
              int j0 = j;
              int j2 = 2*width-j-straddle+1;
              int cc0 = CCTK_GFINDEX3D(cctkGH,i,j0,k);
              int cc2 = CCTK_GFINDEX3D(cctkGH,i,j2,k);
              var[cc0] = sym*var[cc2];
            }
          }
        }
      } else if(face == 3) {
        for(int k=0;k<cctk_lsh[2];k++) {
          for(int j=0;j<width;j++) {
            for(int i=0;i<cctk_lsh[0];i++) {
              int j0 = cctk_lsh[1]-1-j;
              int j2 = cctk_lsh[1]-2*width+j-2+straddle;
              int cc0 = CCTK_GFINDEX3D(cctkGH,i,j0,k);
              int cc2 = CCTK_GFINDEX3D(cctkGH,i,j2,k);
              var[cc0] = sym*var[cc2];
            }
          }
        }
      } else if(face == 4) {
        for(int k=0;k<width;k++) {
          for(int j=0;j<cctk_lsh[1];j++) {
            for(int i=0;i<cctk_lsh[0];i++) {
              int k0 = k;
              int k2 = 2*width-k-straddle+1;
              int cc0 = CCTK_GFINDEX3D(cctkGH,i,j,k0);
              int cc2 = CCTK_GFINDEX3D(cctkGH,i,j,k2);
              var[cc0] = sym*var[cc2];
            }
          }
        }
      } else if(face == 5) {
        for(int k=0;k<width;k++) {
          for(int j=0;j<cctk_lsh[1];j++) {
            for(int i=0;i<cctk_lsh[0];i++) {
              int k0 = cctk_lsh[2]-1-k;
              int k2 = cctk_lsh[2]-2*width+k-2+straddle;
              int cc0 = CCTK_GFINDEX3D(cctkGH,i,j,k0);
              int cc2 = CCTK_GFINDEX3D(cctkGH,i,j,k2);
              var[cc0] = sym*var[cc2];
            }
          }
        }
      } else {
        abort();
      }
    }
  }
  return 0;
}

void presync_registerboundary(CCTK_ARGUMENTS)
{
  Boundary_RegisterPhysicalBC(cctkGH,(boundary_function)fun_bf2<1,1>,"symmetry");
  Boundary_RegisterPhysicalBC(cctkGH,(boundary_function)fun_bf2<-1,1>,"antisymmetry");
  Boundary_RegisterPhysicalBC(cctkGH,(boundary_function)fun_stwave,"zero");
}

void presync_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS

  Boundary_SelectGroupForBC(cctkGH,
    CCTK_ALL_FACES, 1,
   -1 /* no table */, "PresyncWave::evo_vars",
   BCtype);
}

void energy_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS

  Boundary_SelectGroupForBC(cctkGH,
    CCTK_ALL_FACES, 1,
   -1 /* no table */, "PresyncWave::wave_energy",
   BCtype);
}
