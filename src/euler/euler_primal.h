// Euler physics funcitons
#ifndef EULER_PRIMAL_H
#define EULER_PRIMAL_H

// used for defining statically sized arrays
#define MAX_DOFPERNODE 5

#include <cmath>  // sqrt
#include <complex> // complex sqrt
#include "utils.h"
#include "ADTestConfig.h"


namespace Ticon{


//-----------------------------------------------------------------------------
// declarations

template <typename Tsol, typename Tmsh, typename Tres>
void calcSAT(Tsol* roe_vars, Tsol* dq, Tmsh* nrm, Tres* sat);


template <typename Tsol, typename Tmsh, typename Tres>
void calcEulerFlux (Tsol* q, Tres* aux_vars, Tmsh* dir, Tres* F);

template <typename Tsol>
Tsol calcPressure(Tsol* q);



//-----------------------------------------------------------------------------
// definitions of templates functions

#define cdouble Tsol  // datatype for constants, allows easy switching back and
                      // forth

template <typename Tsol, typename Tmsh, typename Tres>
void RoeSolver(Tsol* q, Tsol* qg, Tres* aux_vars, Tmsh* nrm, Tres* flux)
{
  // SAT terms are used for ensuring consistency with the physical problem. Its
  // similar to upwinding which adds dissipation to the problem. SATs on the
  // boundary can be thought of as having two overlapping nodes and because of
  // the discontinuous nature of SBP adds some dissipation.

  const int numDofPerNode = 4;  // TODO: get from Params
  // Declaring constants
  cdouble d1_0 = 1.0;
  cdouble d0_0 = 0.0;
  cdouble d0_5 = 0.5;
  cdouble tau = 1.0;
  cdouble gamma = 1.4;
  cdouble gami = 0.4;
//  gamma = params.gamma
//  gami = params.gamma_1
  cdouble sat_fac = 1;  // multiplier for SAT term

  // Begin main executuion
  Tmsh nx = nrm[0];
  Tmsh ny = nrm[1];

  // Compute the Roe Averaged states
  // The left state of Roe are the actual solution variables
  Tsol fac = d1_0/q[0];
  Tsol uL = q[1]*fac; Tsol vL = q[2]*fac;
  Tsol phi = d0_5*(uL*uL + vL*vL);
  Tsol HL = gamma*q[3]*fac - gami*phi; // Total enthalpy, H = e + 0.5*(u^2 + v^2) + p/rho,
                                 // where e is the internal energy per unit mass

  // The right side of the Roe solver comprises the boundary conditions
       fac = d1_0/qg[0];
  Tsol uR = qg[1]*fac; Tsol vR = qg[2]*fac;
       phi = d0_5*(uR*uR + vR*vR);
  Tsol HR = gamma*qg[3]*fac - gami*phi; // Total Enthalpy

  // Averaged states
  Tsol sqL = sqrt(q[0]);
  Tsol sqR = sqrt(qg[0]);
       fac = d1_0/(sqL + sqR);
  Tsol u = (sqL*uL + sqR*uR)*fac;
  Tsol v = (sqL*vL + sqR*vR)*fac;
  Tsol H = (sqL*HL + sqR*HR)*fac;

//  dq = params.v_vals2 // zeros(Tsol, 4)
  Tsol dq[MAX_DOFPERNODE];
  for (int i=0; i < numDofPerNode; ++i)
    dq[i] = q[i] - qg[i];

//  sat = params.sat_vals
  Tres sat[MAX_DOFPERNODE];
  Tsol roe_vars[3];
//  roe_vars = params.roe_vars
  roe_vars[0] = u;
  roe_vars[1] = v;
  roe_vars[2] = H;
  calcSAT(roe_vars, dq, nrm, sat);

  Tres euler_flux[MAX_DOFPERNODE];
//  euler_flux = params.flux_vals1;
  // calculate Euler flux in wall normal directiona
  // because edge numbering is rather arbitary, any memory access is likely to
  // be a cache miss, so we recalculate the Euler flux
//  v_vals = params.q_vals
//  nrm2 = params.nrm
  Tmsh nrm2[1];
  nrm2[0] = nx;   // why are we assigning to nrm2?
  nrm2[1] = ny;

//  convertFromNaturalToWorkingVars(params, q, v_vals)
  calcEulerFlux(q, aux_vars, nrm2, euler_flux);

  for (int i=0; i < numDofPerNode; ++i)
    flux[i] = (sat_fac*sat[i] + euler_flux[i]);
}


// why is sat Tsol rather than Tres?
template <typename Tsol, typename Tmsh, typename Tres>
void calcSAT(Tsol* roe_vars, Tsol* dq, Tmsh* nrm, Tres* sat)
{
// roe_vars = [u, v, H] at Roe average 

  const int numDofPerNode = 4;

  // SAT parameters
  Tsol sat_Vn = 0.025;
  Tsol sat_Vl = 0.025;
  cdouble tau = 1.0;

  Tsol u = roe_vars[0];
  Tsol v = roe_vars[1];
  Tsol H = roe_vars[2];

  const cdouble gami = 0.4; // TODO: params

  // Begin main executuion
  Tmsh nx = nrm[0];
  Tmsh ny = nrm[1];

  Tmsh dA = sqrt(nx*nx + ny*ny);

  Tres Un = u*nx + v*ny; // Normal Velocity

  Tsol phi = 0.5*(u*u + v*v);

  Tsol a = sqrt(gami*(H - phi)); // speed of sound

  Tres lambda1 = Un + dA*a;
  Tres lambda2 = Un - dA*a;
  Tres lambda3 = Un;

  Tres rhoA = absvalue(Un) + dA*a;

  // Compute Eigen Values of the Flux Jacobian
  // The eigen values calculated above cannot be used directly. Near stagnation
  // points lambda3 approaches zero while near sonic lines lambda1 and lambda2
  // approach zero. This has a possibility of creating numerical difficulties.
  // As a result, the eigen values are limited by the following expressions.

  lambda1 = 0.5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) - lambda1);
  lambda2 = 0.5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) - lambda2);
  lambda3 = 0.5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) - lambda3);


  Tsol dq1 = dq[0];
  Tsol dq2 = dq[1];
  Tsol dq3 = dq[2];
  Tsol dq4 = dq[3];

  sat[0] = lambda3*dq1;
  sat[1] = lambda3*dq2;
  sat[2] = lambda3*dq3;
  sat[3] = lambda3*dq4;


//  E1dq = params.res_vals1
//  E2dq = params.res_vals2
  Tres E1dq[MAX_DOFPERNODE];
  Tres E2dq[MAX_DOFPERNODE];

  //-- get E1*dq
  E1dq[0] = phi*dq1 - u*dq2 - v*dq3 + dq4;
  E1dq[1] = E1dq[0]*u;
  E1dq[2] = E1dq[0]*v;
  E1dq[3] = E1dq[0]*H;

  //-- get E2*dq
  E2dq[0] = 0.0;
  E2dq[1] = -Un*dq1 + nx*dq2 + ny*dq3;
  E2dq[2] = E2dq[1]*ny;
  E2dq[3] = E2dq[1]*Un;
  E2dq[1] = E2dq[1]*nx;

  //-- add to sat
  Tres tmp1 = 0.5*(lambda1 + lambda2) - lambda3;
  Tres tmp2 = gami/(a*a);
  Tres tmp3 = 1.0/(dA*dA);

  for (int i=0; i < numDofPerNode; ++i)
    sat[i] = sat[i] + tmp1*(tmp2*E1dq[i] + tmp3*E2dq[i]);

  //-- get E3*dq
  E1dq[0] = -Un*dq1 + nx*dq2 + ny*dq3;
  E1dq[1] = E1dq[0]*u;
  E1dq[2] = E1dq[0]*v;
  E1dq[3] = E1dq[0]*H;

  //-- get E4*dq
  E2dq[0] = 0.0;
  E2dq[1] = phi*dq1 - u*dq2 - v*dq3 + dq4;
  E2dq[2] = E2dq[1]*ny;
  E2dq[3] = E2dq[1]*Un;
  E2dq[1] = E2dq[1]*nx;

  //-- add to sat
  tmp1 = 0.5*(lambda1 - lambda2)/(dA*a);
  for (int i=0; i < numDofPerNode; ++i)
    sat[i] = sat[i] + tmp1*(E1dq[i] + gami*E2dq[i]);
}

template <typename Tsol, typename Tmsh, typename Tres>
void calcEulerFlux (Tsol* q, Tres* aux_vars, Tmsh* dir, Tres* F)
{
// calculates the Euler flux in a particular direction at a point
// eqn is the equation type
// q is the vector (of length 4), of the conservative variables at the point
// aux_vars is the vector of auxiliary variables at the point
// dir is a vector of length 2 that specifies the direction
// F is populated with the flux (is a vector of length 4)
// 2D  only


  Tsol press = calcPressure(q);
//  press = getPressure(aux_vars)
//  press = @getPressure(aux_vars)
  Tres U = (q[1]*dir[0] + q[2]*dir[1])/q[0];
  F[0] = q[0]*U;
  F[1] = q[1]*U + dir[0]*press;
  F[2] = q[2]*U + dir[1]*press;
  F[3] = (q[3] + press)*U;

}

template <typename Tsol>
Tsol calcPressure(Tsol* q)
{
  // calculate pressure for a node
  // q is a vector of length 4 of the conservative variables

  const cdouble gamma_1 = 0.4;
  return (gamma_1)*(q[3] - 0.5*(q[1]*q[1] + q[2]*q[2])/q[0]);
}

// Instantiate commonly used templates?

} // namespace
#endif
