/*
**
**   rolldown.c
**
**   given a starting point on an equilibrium, roll down the hill and
**   output the points as a function of psi
**
**   2/5/99 -- Darren T. Garnier, Columbia University
**
**
**
**
**
*/

#ifndef _ROLLDOWN_
#define _ROLLDOWN_ 1

#ifdef __cplusplus
extern "C" {
#endif

void RollDownHill(TOKAMAK *td, double r0, double z0,
                  int n, double *psi, double *rr, double *zz);
/* psi, r, and z should be vectors that are already set up going from 0 - n-1 */
/* r0 and z0 are the initial positions */

void GetRandV(TOKAMAK *td, double r0, double z0,
                  int n, double *psi, double *rr, double *zz, double *vv);
/* psi, r, and z should be vectors that are already set up going from 0 - n-1 */
/* r0 and z0 are the initial positions */

void GetRandVfromPsi(TOKAMAK *td, int n,
                     double *psi, double *rr, double *zz, double *vv);
/* psi, r, and z should be vectors that are already set up going from 0 - n-1 */
/* initial position found from psi[0] */

#ifdef __cplusplus
}
#endif

#endif
