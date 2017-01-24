#include <cctk.h>
#include <iostream>

extern "C" void Carpet_SetValidRegion(int,int,int);
extern "C" int Carpet_GetValidRegion(int,int);
extern "C" void ShowValid();

void UpdateValidForMoLInitialCopy(CCTK_ARGUMENTS)
{
  _DECLARE_CCTK_ARGUMENTS;
  const int n  = CCTK_NumVars();
  for(int varindex=0;varindex<n;varindex++) {
    int rhsindex = MoLQueryEvolvedRHS(varindex);
    if(rhsindex >= 0) {
      int mask1 = Carpet_GetValidRegion(varindex,1);
      int mask2 = Carpet_GetValidRegion(rhsindex,1);
      Carpet_SetValidRegion(varindex,0,mask1);
      Carpet_SetValidRegion(rhsindex,0,mask2);
    }
  }
}

void UpdateValidForMoLAdd(CCTK_ARGUMENTS)
{
  _DECLARE_CCTK_ARGUMENTS;
  const int n  = CCTK_NumVars();
  for(int varindex=0;varindex<n;varindex++) {
    int rhsindex = MoLQueryEvolvedRHS(varindex);
    if(rhsindex >= 0) {
      int mask1 = Carpet_GetValidRegion(varindex,0);
      int mask2 = Carpet_GetValidRegion(varindex,1);
      int mask3 = Carpet_GetValidRegion(rhsindex,0);
      int mask = mask1 & mask2 & mask3;
      Carpet_SetValidRegion(varindex,0,mask);
    }
  }
}
