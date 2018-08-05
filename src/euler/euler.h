// main header for all Euler functions

#ifndef EULER_H
#define EULER_H

template <typename Tsol, typename Tmsh, typename Tres>
void RoeSolver(Tsol* q, Tsol* qg, Tres* aux_vars, Tmsh* nrm, Tres* flux);


#include "euler_primal.h"
#include "euler_diff.h"

#ifdef ENABLE_CODIPACK
  #include "euler_codipack.h"
#endif

#ifdef ENABLE_ADEPT
  #include "euler_adept.h"
#endif

#endif
