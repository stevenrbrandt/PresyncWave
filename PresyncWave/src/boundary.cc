#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <iostream>

typedef void (*boundary_function)(
  const cGH *cctkGH,
  int num_vars,
  int *var_indices,
  int *faces,
  int *widths,
  int *table_handles);

// From PreSync
extern "C" void Carpet_RegisterPhysicalBC(
    const cGH *cctkGH,
    boundary_function func,
    const char *bc_name,
    int before);

extern "C"
void Carpet_SelectGroupForBC(
    const cGH *cctkGH,
    int faces,
    int width,
    int table_handle,
    const char *group_name,
    const char *bc_name);

extern "C" void Boundary_RegisterPhysicalBC(CCTK_ARGUMENTS,boundary_function,const char *);

// 1 = symmetric, -1 = antisymmetric
// 1 = straddle: point is a grid point
// 0 = straddle: point is between grid points.

template<int sym,int straddle>
void fun_blank(
  const cGH *cctkGH,
  int num_vars,
  int *var_indices,
  int *faces,
  int *widths,
  int *table_handles) {
}

template<int sym,int straddle>
void fun_bf2(
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
    std::cout << "BOUNDARY: " << CCTK_FullName(var_index) << std::endl;
    for(int face=0;face<6;face++) {

      if(!cctk_bbox[face])
        continue; // we aren't on a physical boundary

      int width = widths[vi];
      assert(width > 0);
      CCTK_REAL *var = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH,0,var_index);
      if(face == 0) {
        for(int k=0;k<cctk_lsh[2];k++) {
          for(int j=0;j<cctk_lsh[1];j++) {
            for(int i=0;i<width;i++) {
              int i0 = i;
              int i2 = 2*width-i+1-straddle;
              if(j==0 && k==0) {
                //std::cout << "i0=" << i0 << " i2=" << i2 << std::endl;
              }
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
              if(j==0 && k==0) {
                //std::cout << "i0=" << i0 << " i2=" << i2 << std::endl;
              }
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
}

void fun_bf(
  const cGH *cctkGH,
  int num_vars,
  int *var_indices,
  int *faces,
  int *widths,
  int *table_handles) {
  _DECLARE_CCTK_ARGUMENTS;
  
  std::cout << "fun: num_vars=" << num_vars << std::endl;
  for(int vi=0;vi<num_vars;vi++) {
    assert(widths[vi] > 0);
    int var_index = var_indices[vi];
    for(int face=0;face<6;face++) {

      if(!cctk_bbox[face])
        continue; // we aren't on a physical boundary

      int dim = face/2; // 0,1=>0 2,3=>1, 4,5=>2
      int dir[3] = {0,0,0};
      if(face % 2 == 1) {
        dir[dim] = -1;
      } else {
        dir[dim] = 1;
      }

      int zero = CCTK_GFINDEX3D(cctkGH,0,0,0);
      int max_index = CCTK_GFINDEX3D(cctkGH,cctk_lsh[0],cctk_lsh[1],cctk_lsh[2]);
      int lo[3], hi[3], del[3], one[3];
      for(int d=0;d<3;d++) {
        if(d == dim) {
          if(face % 2 == 1) {
            lo[d] = cctk_lsh[d]-widths[vi];
            hi[d] = cctk_lsh[d];
          } else {
            lo[d] = 0;
            hi[d] = widths[vi];
          }
        } else {
          lo[d] = 0;
          hi[d] = cctk_lsh[d];
        }
        // Compute a unit vector in the d dimension
        for(int k=0;k<3;k++)
          one[k] = k==d ? 1 : 0;
        del[d] = CCTK_GFINDEX3D(cctkGH,one[0],one[1],one[2]) - zero;
      }

      // This ought to be true
      assert(del[0] == 1);

      CCTK_REAL *var = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH,0,var_index);
      //assert(var != 0);
      if(var == 0) {
        std::cout << "NOT updating boundary for " << CCTK_FullName(var_index) << "(" << var_index << ") NULL PTR!" << std::endl;
        abort();
        return;
      } else {
        std::cout << "updating boundary for " << CCTK_FullName(var_index) << "(" << var_index << "), face=" << face << std::endl;
      }

      int ind[3];
      for(ind[2]=lo[2];ind[2] < hi[2];ind[2]++) {
        for(ind[1]=lo[1];ind[1] < hi[1];ind[1]++) {
          int cc0 = zero + del[1]*ind[1] + del[2]*ind[2];
          for(ind[0]=lo[0];ind[0] < hi[0];ind[0]++) {
            int cc = cc0 + ind[0];

            int ci = cc;
            if(face % 2 == 1) {
              ci += 2*(lo[dim]-ind[dim])*del[dim]; 
            } else {
              ci += 2*(hi[dim]-ind[dim])*del[dim];
              assert(lo[dim]==0);
            }

            assert(cc < max_index);
            if(ci >= max_index || ci < zero) {
              std::cout << "dim=" << dim << " w=" << widths[vi] << " ind=" << ind[dim] << " lo=" << lo[dim] << " hi=" << hi[dim] << " up=" << cctk_lsh[dim] << std::endl;
            }
            assert(ci < max_index);
            assert(cc >= zero);
            assert(ci >= zero);
            std::cout << " face=" << face;
            std::cout << " var(" << ind[0] << "," << ind[1] << "," << ind[2] << ")";
            std::cout << " var2(";
            for(int i=0;i<3;i++) {
              if(i != 0)   std::cout << ",";
              if(dim == i) std::cout << 2*lo[dim]-ind[dim];
              else         std::cout << ind[i];
            }
            std::cout << ")" << std::endl;
            var[cc] = var[ci]; // symmetry boundary condition
          }
        }
      }
    }
  }
}

extern "C"
void presync_registerboundary(CCTK_ARGUMENTS)
{
  _DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  CCTK_INT group, rhs;

  std::cout << "Register Boundary Conditions" << std::endl;

  Carpet_RegisterPhysicalBC(cctkGH,fun_bf2<1,1>,"symmetry",1);
  Carpet_RegisterPhysicalBC(cctkGH,fun_bf2<-1,1>,"antisymmetry",1);
//  Carpet_RegisterPhysicalBC(cctkGH,fun_zero,"antisymmetry",1);
  int w = 1;

  Carpet_SelectGroupForBC(cctkGH,
    CCTK_ALL_FACES, w,
   -1 /* no table */, "PresyncWave::phi_g",
   "antisymmetry");

  Carpet_SelectGroupForBC(cctkGH,
    CCTK_ALL_FACES, w,
   -1 /* no table */, "PresyncWave::phi_rhs_g",
   "antisymmetry");

  Carpet_SelectGroupForBC(cctkGH,
    CCTK_ALL_FACES, w,
   -1 /* no table */, "PresyncWave::psi_g",
   "antisymmetry");

  Carpet_SelectGroupForBC(cctkGH,
    CCTK_ALL_FACES, w,
   -1 /* no table */, "PresyncWave::psi_rhs_g",
   "antisymmetry");
}
