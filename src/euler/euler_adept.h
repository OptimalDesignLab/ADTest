// using Adept for differentiation

#ifndef EULER_ADEPT_H
#define EULER_ADEPT_H

#ifdef ENABLE_ADEPT

#include "adept.h"
#include "euler_primal.h"

namespace Ticon {

template <typename Tsol, typename Tmsh, typename Tres>
void RoeSolver_diff_adept( adept::Stack& stack, Tsol* q, Tsol* qg, Tmsh* nrm,
                           Tres* fluxL_dot, Tres* fluxR_dot)
{
  // Tsol = Tres should be double, because this function internally handles
  // conversion to the AD datatype
  using adept::adouble;
//  static adept::Stack stack;
  adouble flux[4];
  adouble _q[4] = {q[0], q[1], q[2], q[3]};
  adouble _qg[4] = {qg[0], qg[1], qg[2], qg[3]};
  adouble aux_vars[0];

  stack.new_recording();
  RoeSolver(_q, _qg, aux_vars, nrm, flux);

  stack.independent(_q, 4);
  stack.independent(_qg, 4);
  stack.dependent(flux, 4);

  double jac_tmp[4*4*2]; // TODO: avoid allocating this temporary
  stack.jacobian(jac_tmp);

  for (int i=0; i < 16; ++i)
  {
    fluxL_dot[i] = jac_tmp[i];
    fluxR_dot[i] = jac_tmp[i+16];
  }
}

}

#endif // end if ENABLE_ADEPT
#endif  // end header guard
