/*
 *
 * Date: 
 *   Revised: May 08, 2014
 * 
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "densitycut.h"

/* Register all the functions */
static const R_CallMethodDef callMethods[] = {
  {"sexp_log_sum_exp", (DL_FUNC) &sexp_log_sum_exp, 1},
  {NULL, NULL, 0}
};

void R_init_densitycut(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
