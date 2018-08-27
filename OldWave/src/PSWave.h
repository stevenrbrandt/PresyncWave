#ifndef _PSWave_H_
#define _PSWave_H_

typedef CCTK_INT (*phys_bc_fn_ptr)(CCTK_POINTER_TO_CONST, const CCTK_INT,
                                   const CCTK_INT *, const CCTK_INT *,
                                   const CCTK_INT *, const CCTK_INT *);

#ifdef __cplusplus
extern "C" {
#endif

/* prototype for routine registered as providing 'zero' boundary condition */
CCTK_INT fun_stwave(const cGH *cctkGH, CCTK_INT num_vars, CCTK_INT *var_indices,
                  CCTK_INT *faces, CCTK_INT *widths, CCTK_INT *table_handles);

/* prototype for routine registered as providing 'symmetry' and 'anti-symmetry' boundary conditions */
CCTK_INT fun_bf2(const cGH *cctkGH, CCTK_INT num_vars, CCTK_INT *var_indices,
                  CCTK_INT *faces, CCTK_INT *widths, CCTK_INT *table_handles);

#ifdef __cplusplus
}
#endif

#endif /* _PSWave_H_ */
