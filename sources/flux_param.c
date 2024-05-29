// flux_param.c
// this is the interface file to IDL.

#include <stdio.h>
#include "export.h"

#include "psigrid.h"
#include "plasma.h"
#include "tokamak.h"
#include "separatrix.h"
#include "PlasmaBoundary.h"
#include "GetPlasmaParameters.h"
#include "GetFluxMoments.h"
#include "GetFluxParameters.h"


#define MIN_FP_POS_ARGS 3
#define MAX_FP_POS_ARGS 8

__declspec(export) int IDL_Load(void);
void flux_param_idl(int argc, IDL_VPTR argv[], char *argk);

#define ARRLEN(arg) sizeof(arg)/sizeof(arg[0])

int IDL_Load(void)
{
  /* These tables contain information on the functions and procedures
   * that make up the TESTMODULE DLM. The information contained in these
   * tables must be identical to that contained in testmodule.dlm.
   */
//  static IDL_SYSFUN_DEF function_addr[] = {
//    { testfun, "TESTFUN", 0, IDL_MAXPARAMS, 0},
//  };

  static IDL_SYSFUN_DEF procedure_addr[] = {
    { (IDL_FUN_RET) flux_param_idl, "FLUX_PARAM",
        MIN_FP_POS_ARGS, MAX_FP_POS_ARGS, IDL_SYSFUN_DEF_F_KEYWORDS},
  };

  /* Create a message block to hold our messages. Save its handle where
   * the other routines can access it. */
//  if (!(msg_block = IDL_MessageDefineBlock("Testmodule", ARRLEN(msg_arr),
//					           msg_arr))) return IDL_FALSE;

  /* Register our routine. The routines must be specified exactly the same
   * as in testmodule.dlm. */
  return // IDL_AddSystemRoutine(function_addr, TRUE, ARRLEN(function_addr)) &&
         IDL_AddSystemRoutine(procedure_addr, FALSE, ARRLEN(procedure_addr)) &&
         1;
}


IDL_VPTR my_function(int kargc, IDL_VPTR kargv[], char *argk);


void flux_param_idl(int kargc, IDL_VPTR kargv[], char *argk)
{
   IDL_VPTR pg_v, rg_v, zg_v, temp_v;
   double *p, *r, *z, *t;
   TOKAMAK *td;
   PSIGRID *pg;
   PLASMA *pl;
   SEPARATRIX *s;
   LIMITER *l;
   int i,j,nsize;
   IDL_LONG dim[3];

   int argc;
   IDL_VPTR argv[MAX_FP_POS_ARGS];

   static int w_psiaxis, w_psilcfs;
   static IDL_VPTR v_psiaxis, v_psilcfs;
   static IDL_LONG npsi;
   static double psimax;

   static IDL_KW_PAR kw_pars[] = { IDL_KW_FAST_SCAN,
      { "NPSI",     IDL_TYP_LONG,   1, IDL_KW_ZERO,         NULL,    IDL_CHARA(npsi)      },
      { "PSIAXIS",  IDL_TYP_DOUBLE, 1, IDL_KW_OUT,     &w_psiaxis ,  IDL_CHARA(v_psiaxis) },
      { "PSILCFS",  IDL_TYP_DOUBLE, 1, IDL_KW_OUT,     &w_psilcfs ,  IDL_CHARA(v_psilcfs) },
      { "PSIMAX",   IDL_TYP_DOUBLE, 1, IDL_KW_ZERO,         NULL,    IDL_CHARA(psimax)    },
      { NULL }
   };


   IDL_KWCleanup(IDL_KW_MARK);

   argc = IDL_KWGetParams(kargc, kargv, argk, kw_pars, argv, 1);


   /* flux param have at least 3 inputs  */

   if (argc <= 3)
      IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,"Too few arguments.");


   pg_v = argv[0];
   rg_v = argv[1];
   zg_v = argv[2];

   IDL_ENSURE_ARRAY(pg_v);
   IDL_ENSURE_SIMPLE(pg_v);
   if ((pg_v->value.arr->n_dim  != 2 )  ||
       (pg_v->value.arr->dim[0] <= 4 )  ||
       (pg_v->value.arr->dim[0] <= 4 )  ||
       (pg_v->value.arr->dim[0] != pg_v->value.arr->dim[1] ))
      IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,"PSI must be square matrix.");

   IDL_ENSURE_ARRAY(rg_v);
   IDL_ENSURE_SIMPLE(rg_v);
   IDL_ENSURE_ARRAY(zg_v);
   IDL_ENSURE_SIMPLE(zg_v);
   if ((rg_v->value.arr->n_elts != 3 ) ||
       (zg_v->value.arr->n_elts != 3 ) )
      IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_LONGJMP,"R and Z are 3 element vectors. (min,max,axis)");

   if (pg_v->type != IDL_TYP_DOUBLE) pg_v = IDL_CvtDbl(1,&pg_v);
   if (rg_v->type != IDL_TYP_DOUBLE) rg_v = IDL_CvtDbl(1,&rg_v);
   if (zg_v->type != IDL_TYP_DOUBLE) zg_v = IDL_CvtDbl(1,&zg_v);

   /* now we have psi as a double array and vectors of rgrid and zgrid */

   p = (double *) pg_v->value.arr->data;
   r = (double *) rg_v->value.arr->data;
   z = (double *) zg_v->value.arr->data;

   td = new_Tokamak();
   pg = td->PsiGrid;
   pl = td->Plasma;

   /* pre-init default fixing. */
   nsize = pg_v->value.arr->dim[0];

   pg->Nsize = nsize-1;
   pg->Symmetric = 0;
   pg->Xmin = r[0];
   pg->Xmax = r[1];
   pg->Zmin = z[0];
   pg->Zmax = z[1];


   pl->Nsize = nsize;

   pl->a0 = 3.0;
   pl->Ip0 = 0;         /* no plasma current. */
   pl->B0 = 0;          /* no toroidal field */
   if (npsi < 1) npsi = 25;
   pl->NumPsiPts = npsi;
   if (psimax < .1) psimax = .98;
   pl->PsiXmax = psimax;  /* get a little closer to the edge */
   pl->XMagAxis = pl->R0 = r[2];
   pl->ZMagAxis = pl->Z0 = z[2];

   td->NumLimiters= 2;
   td->NumSeps=1;

   init_Tokamak(td);

   td->Seps[0] = s = new_Separatrix();
   s->X1 = .5;
   s->Z1 = .9;
   s->X2 = 1.7;
   s->Z2 = 1.5;
   s->Enabled = 1;


   td->Limiters[0] = l = new_Limiter();
   l->X1 = 0.2438;   /* this is just a line from inside to outside major radius */
   l->Z1 = 0.0;      /* reference LDX-MIT-AZ-032398-01 */
   l->X2 = 0.5642;
   l->Z2 = 0.0;
   l->Enabled = 2;   /* special DIPOLE limiter */

   td->Limiters[1] = l = new_Limiter();
   l->X1 = 0.5642;   /* this is just a line from inside to outside major radius */
   l->Z1 = 0.0;      /* reference LDX-MIT-AZ-032398-01 */
   l->X2 = 0.2438;
   l->Z2 = 0.0;
   l->Enabled = 2;   /* special DIPOLE limiter */

   /* now lets copy our data */

   for (i=0; i<nsize; i++)
       for (j=0; j<nsize; j++)
          pg->Psi[j][i] = p[i*nsize+j];  /* do idl row/colum reverse ? */



   PlasmaBoundary(td);

   /* compute the gradPsi parameters so we can do the integrations */
   /* from GetPlasmaParameters */

   //GetGradPsi(td);
   //GetFluxParameters(td);

   /* here we should really find the two max and min flux for matching.... */
   /* we should assume that min psi is near the coil */
   /* and do the checks against limiters */
   /* so lets assume we have a boundary in R and Z for limiters inner and outer..... */


   GetPlasmaParameters(td);


   /* now we shall return all parameters....
    * lets return the following....
    * array of params = PsiXmax, PsiAxis, PsiLim, XMagAxis, ZMagAxis
    * psi array
    * dv/dpsi(psi)
    * v(psi)
    * B^2(psi)
    * Well(psi)
    */

   /* first array of parameters... */

   dim[0] = pl->NumPsiPts;
   t = (double *)IDL_MakeTempArray(IDL_TYP_DOUBLE, 1, dim, IDL_BARR_INI_NOP, &temp_v);
   for (i=0; i < dim[0]; i++)
        t[i] = pg->PsiAxis + i * pl->PsiXmax/(pl->NumPsiPts-1.0) * pg->DelPsi;

   IDL_VarCopy(temp_v,argv[3]);

   /* now the profiles.... */

   if (argc <= 4) {
       free_Tokamak(td, TRUE);
       return;
   }

   /* d V / d P s i */

   dim[0] = pl->NumPsiPts;
   t = (double *)IDL_MakeTempArray(IDL_TYP_DOUBLE,1,dim,IDL_BARR_INI_NOP, &temp_v);
   for (i=0; i < dim[0]; i++) t[i]=pl->Volp_pr[i];
   IDL_VarCopy(temp_v,argv[4]);

   if (argc <= 5) {
       free_Tokamak(td, TRUE);
       return;
   }

   /* V o l u m e */

   dim[0] = pl->NumPsiPts;
   t = (double *)IDL_MakeTempArray(IDL_TYP_DOUBLE,1,dim,IDL_BARR_INI_NOP, &temp_v);
   for (i=0; i < dim[0]; i++) t[i]=pl->Vol_pr[i];
   IDL_VarCopy(temp_v,argv[5]);

   if (argc <= 6) {
       free_Tokamak(td, TRUE);
       return;
   }

   /* B ^ 2 */

   dim[0] = pl->NumPsiPts;
   t = (double *)IDL_MakeTempArray(IDL_TYP_DOUBLE,1,dim,IDL_BARR_INI_NOP, &temp_v);
   for (i=0; i < dim[0]; i++) t[i]=pl->B2_pr[i];
   IDL_VarCopy(temp_v,argv[6]);

   if (argc <= 7) {
       free_Tokamak(td, TRUE);
       return;
   }

   /* M a g n e t i c   W e l l */

   dim[0] = pl->NumPsiPts;
   t = (double *)IDL_MakeTempArray(IDL_TYP_DOUBLE,1,dim,IDL_BARR_INI_NOP, &temp_v);
   for (i=0; i < dim[0]; i++) t[i]=pl->Well_pr[i];
   IDL_VarCopy(temp_v,argv[7]);


   if (w_psiaxis != 0)
      IDL_StoreScalar(v_psiaxis,IDL_TYP_DOUBLE, (IDL_ALLTYPES *) &pg->PsiAxis);
   if (w_psilcfs != 0)
      IDL_StoreScalar(v_psilcfs,IDL_TYP_DOUBLE, (IDL_ALLTYPES *) &pg->PsiLim);




   free_Tokamak(td, TRUE);

   IDL_KWCleanup(IDL_KW_CLEAN);

   return;

}
