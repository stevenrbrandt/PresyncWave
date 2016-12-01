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
void Carpet_SelectVarForBCI(
    const cGH *cctkGH,
    int faces,
    int width,
    int table_handle,
    int var_index,
    const char *bc_name);

extern "C"
void Carpet_ClearBCForVarI(
    const cGH *cctkGH,
    int var_index);

extern "C"
void Carpet_ApplyPhysicalBCS(const cGH *cctkGH);

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
  
  std::cout << "fun: num_vars=" << num_vars << std::endl;
  for(int i=0;i<num_vars;i++) {
    std::cout << "fun:";
    std::cout << " vi=" << CCTK_FullName(var_indices[i]);
    std::cout << " faces=" << faces[i];
    std::cout << " widths=" << widths[i];
    std::cout << " table_handles=" << table_handles[i];
    std::cout << std::endl;
  }
}
void fun_bf2(
  const cGH *cctkGH,
  int num_vars,
  int *var_indices,
  int *faces,
  int *widths,
  int *table_handles) {
  
  std::cout << "fun2: num_vars=" << num_vars << std::endl;
  for(int i=0;i<num_vars;i++) {
    std::cout << "fun2:";
    std::cout << " vi=" << CCTK_FullName(var_indices[i]);
    std::cout << " faces=" << faces[i];
    std::cout << " widths=" << widths[i];
    std::cout << " table_handles=" << table_handles[i];
    std::cout << std::endl;
  }
}

extern "C"
void presync_registervars (CCTK_ARGUMENTS)
{
  _DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  CCTK_INT ierr = 0, group, rhs;

  std::cout << "Register" << std::endl;

  rhs = CCTK_GroupIndex ("PresyncWave::rhs_vars");
  group = CCTK_GroupIndex ("PresyncWave::evo_vars");

  if (CCTK_IsFunctionAliased ("MoLRegisterEvolvedGroup"))
  {
    ierr |= MoLRegisterEvolvedGroup (group, rhs);
  }
  else
  {
    CCTK_WARN (0, "MoL function not aliased !");
    ierr |= 1;
  }

  Boundary_RegisterPhysicalBC(cctkGH,fun_bf,"fun");
  Carpet_RegisterPhysicalBC(cctkGH,fun_bf,"fun",1);

  Boundary_RegisterPhysicalBC(cctkGH,fun_bf2,"fun2");
  Carpet_RegisterPhysicalBC(cctkGH,fun_bf2,"fun2",1);

  if (CCTK_IsFunctionAliased ("Boundary_SelectGroupForBC")) {
    ierr |= Boundary_SelectGroupForBC(cctkGH,
      8, 1,
      -1 /* no table */, "PresyncWave::rhs_vars",
      "fun");
    Carpet_SelectGroupForBC(cctkGH,
      8, 1,
      -1 /* no table */, "PresyncWave::rhs_vars",
      "fun2");

    ierr |= Boundary_SelectGroupForBC(cctkGH,
      3, 1,
      -1 /* no table */, "PresyncWave::rhs_vars",
      "fun");
    Carpet_SelectGroupForBC(cctkGH,
      3, 1,
      -1 /* no table */, "PresyncWave::rhs_vars",
      "fun2");

    ierr |= Boundary_SelectGroupForBC(cctkGH,
      CCTK_ALL_FACES, 2,
      -1 /* no table */, "PresyncWave::evo_vars",
      "fun");
    Carpet_SelectGroupForBC(cctkGH,
      CCTK_ALL_FACES, 2,
      -1 /* no table */, "PresyncWave::evo_vars",
      "fun2");
  }
  else
  {
    CCTK_WARN (0, "BC function not aliased !");
    ierr |= 1;
  }

  if (ierr)
    CCTK_WARN (0, "Problems registering with MoL or BCs !");

  return;
}
