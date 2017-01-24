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
  CCTK_INT ierr = 0, group, rhs;

  std::cout << "Register Boundary Conditions" << std::endl;

  Carpet_RegisterPhysicalBC(cctkGH,fun_bf,"symmetry",1);

  Carpet_SelectGroupForBC(cctkGH,
    CCTK_ALL_FACES, 1,
   -1 /* no table */, "PresyncWave::rhs_vars",
   "symmetry");

  Carpet_SelectGroupForBC(cctkGH,
    CCTK_ALL_FACES, 1,
   -1 /* no table */, "PresyncWave::evo_vars",
   "symmetry");
}