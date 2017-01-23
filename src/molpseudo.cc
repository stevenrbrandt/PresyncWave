#include <cctk.h>
#include <iostream>

extern "C" void SetValidRegion(int,int,int);
extern "C" int GetValidRegion(int,int);
extern "C" void ShowValid();

void UpdateValidForMoLInitialCopy(CCTK_ARGUMENTS)
{
  _DECLARE_CCTK_ARGUMENTS;
  const int n  = CCTK_NumVars();
  for(int varindex=0;varindex<n;varindex++) {
    int rhsindex = MoLQueryEvolvedRHS(varindex);
    if(rhsindex >= 0) {
      int mask1 = GetValidRegion(varindex,1);
      int mask2 = GetValidRegion(rhsindex,1);
      SetValidRegion(varindex,0,mask1);
      SetValidRegion(rhsindex,0,mask2);
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
      int mask1 = GetValidRegion(varindex,0);
      int mask2 = GetValidRegion(varindex,1);
      int mask3 = GetValidRegion(rhsindex,0);
      int mask = mask1 & mask2 & mask3;
      SetValidRegion(varindex,0,mask);
    }
  }
}
