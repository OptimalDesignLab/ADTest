
#include "codi.hpp"
#include "euler_primal.h"

#ifndef EULER_CODIPACK
#define EULER_CODIPACK

#ifdef ENABLE_CODIPACK

namespace Ticon {

//TODO: use some kind of partial template specialization
template <typename Tsol, typename Tmsh, typename Tres>
void RoeSolver_diff_codi( Tsol* q, Tsol* qg, Tres* aux_vars, Tmsh* nrm,
                          Tres* fluxL_dot, Tres* fluxR_dot)
{
  // Tsol = Tres = codi::RealForwardVec<8> required
  Tres flux[4];
  for (int i=0; i < 4; ++i)
  {
    q[i]. gradient()[i]   = 1.0;
    qg[i].gradient()[i+4] = 1.0;
    flux[i] = 0.0;
  }


  RoeSolver(q, qg, aux_vars, nrm, flux);

  for (int i=0; i < 4; ++i)
    for (int j=0; j < 4; ++j)
    {
      fluxL_dot[ind(i, j, 4, 4)] = flux[i].gradient()[j];
      fluxR_dot[ind(i, j, 4, 4)] = flux[i].gradient()[j+4];
    }

} // function RoeSolver_diff_codi

}  // namespace Ticon

#endif // end if ENABLE_CODIPACK
#endif // end header guard
