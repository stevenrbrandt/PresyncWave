#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Arguments_Checked.h"
#include "cctk_Parameters.h"
#include <iostream>

void oldsync_registervars (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_oldsync_registervars
  DECLARE_CCTK_PARAMETERS
  CCTK_INT ierr = 0, group, rhs;

  std::cout << "Register" << std::endl;

  rhs = CCTK_GroupIndex ("OldWave::rhs_vars");
  group = CCTK_GroupIndex ("OldWave::evo_vars");

  if (CCTK_IsFunctionAliased ("MoLRegisterEvolvedGroup"))
  {
    ierr |= MoLRegisterEvolvedGroup (group, rhs);
  }
  else
  {
    CCTK_WARN (0, "MoL function not aliased !");
    ierr |= 1;
  }

  if (ierr)
    CCTK_WARN (0, "Problems registering with MoL or BCs !");

  return;
}
