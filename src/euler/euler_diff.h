// differentiated Euler functions

#ifndef EULER_DIFF_H
#define EULER_DIFF_H

#include <cmath>  // sqrt
#include <complex> // complex sqrt
#include "utils.h"
#include "ADTestConfig.h"
#include "euler_primal.h"

namespace Ticon{

template <typename Tsol, typename Tmsh, typename Tres>
void RoeSolver_diff( Tsol* q, stol* qg, Tres* aux_vars, Tmsh* nrm,
                     Tres* fluxL_dot, Tres* fluxR_dot)
{
  // SAT terms are used for ensuring consistency with the physical problem. Its
  // similar to upwinding which adds dissipation to the problem. SATs on the
  // boundary can be thought of as having two overlapping nodes and because of
  // the discontinuous nature of SBP adds some dissipation.

  const int numDofPerNode = 4;

  // Declaring constants
  double d1_0 = 1.0;
  double d0_0 = 0.0;
  double d0_5 = 0.5;
  double tau = 1.0;
  const double gamma = 1.4;
  const double gami = 0.4;
//  gamma = params.gamma
//  gami = params.gamma_1
  double sat_fac = 1;  // multiplier for SAT term

  // Begin main executuion
  auto nx = nrm[1];
  auto ny = nrm[2];

  // Compute the Roe Averaged states
  // The left state of Roe are the actual solution variables
  // All the _dot variables are wrt q
  auto fac = d1_0/q[1];
  auto fac_dot1 = -fac*fac;
//  fac_dot1 = -d1_0/(q[1]*q[1])

  auto uL = q[2]*fac;
  auto uL_dot1 = q[2]*fac_dot1;
  auto uL_dot2 = fac;

  auto vL = q[3]*fac;
  auto vL_dot1 = q[3]*fac_dot1;
  auto vL_dot3 = fac;

  auto phi = d0_5*(uL*uL + vL*vL);
  auto phi_dot1 = uL*uL_dot1 + vL*vL_dot1;
  auto phi_dot2 = uL*uL_dot2;
  auto phi_dot3 = vL*vL_dot3;

  auto HL = gamma*q[4]*fac - gami*phi; // Total enthalpy, H = e + 0.5*(u^2 + v^2) + p/rho,
                                 // where e is the internal energy per unit mass
  auto HL_dot1 = gamma*q[4]*fac_dot1 - gami*phi_dot1;
  auto HL_dot2 = -gami*phi_dot2;
  auto HL_dot3 = -gami*phi_dot3;
  auto HL_dot4 = gamma*fac  // q[4]_dot = 1;

  // The right side of the Roe solver comprises the boundary conditions
  // all the _dot variables are wrt qg now
  auto fac = d1_0/qg[1];
  auto fac_dot1 = -fac*fac;
//  fac_dot1 = -d1_0/(qg[1]*qg[1])

  auto uR = qg[2]*fac;
  auto uR_dot1 = qg[2]*fac_dot1;
  auto uR_dot2 = fac;

  auto vR = qg[3]*fac;
  auto vR_dot1 = qg[3]*fac_dot1;
  auto vR_dot3 = fac;

  auto phi = d0_5*(uR*uR + vR*vR);
  auto phi_dot1 = uR*uR_dot1 + vR*vR_dot1;
  auto phi_dot2 = uR*uR_dot2;
  auto phi_dot3 = vR*vR_dot3;


  auto HR = gamma*qg[4]*fac - gami*phi // Total Enthalpy;
  auto HR_dot1 = gamma*qg[4]*fac_dot1 - gami*phi_dot1;
  auto HR_dot2 = -gami*phi_dot2;
  auto HR_dot3 = -gami*phi_dot3;
  auto HR_dot4 = gamma*fac;

  // Averaged states
  auto sqL = sqrt(q[1]);
  auto sqL_dot1 = 0.5/sqL;

  auto sqR = sqrt(qg[1]);
  auto sqR_dot1 = 0.5/sqR;

  fac = d1_0/(sqL + sqR);
  auto t1 = -1/((sqL + sqR)*(sqL + sqR));
  auto fac_dotL1 = t1*sqL_dot1;
  auto fac_dotR1 = t1*sqR_dot1;

  auto u = (sqL*uL + sqR*uR)*fac;
  auto t2 = sqR*uR;
  auto t3 = sqL*uL;
  auto t4 = sqL*fac;
  auto t5 = sqR*fac;
  auto u_dotL1 = t3*fac_dotL1 + sqL*fac*uL_dot1 + uL*fac*sqL_dot1 + t2*fac_dotL1;
  auto u_dotR1 = t2*fac_dotR1 + sqR*fac*uR_dot1 + uR*fac*sqR_dot1 + t3*fac_dotR1;


  auto u_dotL2 = t4*uL_dot2;
  auto u_dotR2 = t5*uR_dot2;

  auto v = (sqL*vL + sqR*vR)*fac;
  auto t2 = sqL*vL;
  auto t3 = sqR*vR;
  auto v_dotL1 = t2*fac_dotL1 + sqL*fac*vL_dot1 + vL*fac*sqL_dot1 + t3*fac_dotL1;
  auto v_dotR1 = t3*fac_dotR1 + sqR*fac*vR_dot1 + vR*fac*sqR_dot1 + t2*fac_dotR1;

  auto v_dotL3 = t4*vL_dot3;
  auto v_dotR3 = t5*vR_dot3;

  auto H = (sqL*HL + sqR*HR)*fac;
  auto t2 = sqL*HL;
  auto t3 = sqR*HR;
  auto H_dotL1 = t2*fac_dotL1 + sqL*fac*HL_dot1 + HL*fac*sqL_dot1 + t3*fac_dotL1;
  auto H_dotR1 = t3*fac_dotR1 + sqR*fac*HR_dot1 + HR*fac*sqR_dot1 + t2*fac_dotR1;
 
  auto H_dotL2 = t4*HL_dot2 ;
  auto H_dotR2 = t5*HR_dot2;

  auto H_dotL3 = t4*HL_dot3;
  auto H_dotR3 = t5*HR_dot3;

  auto H_dotL4 = t4*HL_dot4;
  auto H_dotR4 = t5*HR_dot4;


//  dq = params.v_vals2 // zeros(Tsol, 4)
  Tsol dq[MAX_DOFPERNODE];

  for (int i=0; i < numDofPerNode; ++i)
    dq[i] = q[i] - qg[i];

  // dq_dotL* = 1, dq_dotR* = -1, so omit them

//  roe_vars = params.roe_vars
  Tsol roe_vars [MAX_DOFPERNODE];
  roe_vars[1] = u;
  roe_vars[2] = v;
  roe_vars[3] = H;

//  roe_vars_dot = params.roe_vars_dot
  Tsol roe_vars_dot[16];
  roe_vars_dot[1]  = u_dotL1;
  roe_vars_dot[2]  = u_dotR1;
  roe_vars_dot[3]  = u_dotL2;
  roe_vars_dot[4]  = u_dotR2;

  roe_vars_dot[5]  = v_dotL1;
  roe_vars_dot[6]  = v_dotR1;
  roe_vars_dot[7]  = v_dotL3;
  roe_vars_dot[8]  = v_dotR3;

  roe_vars_dot[9]  = H_dotL1;
  roe_vars_dot[10] = H_dotR1;
  roe_vars_dot[11] = H_dotL2;
  roe_vars_dot[12] = H_dotR2;
  roe_vars_dot[13] = H_dotL3;
  roe_vars_dot[14] = H_dotR3;
  roe_vars_dot[15] = H_dotL4;
  roe_vars_dot[16] = H_dotR4;
  

//  sat = params.sat_vals
  Tres sat[MAX_NUMDOFPERNODE];
//  sat_jacL = params.sat_jacL
//  sat_jacR = params.sat_jacR
  for (int i=0; i < numDofPerNode*numDofPerNode; ++i)
  {
    fluxL_dot[i] = 0;
    fluxR_dot[i] = 0;
  }
  // pass in fluxL_dot and fluxR_dot here, then add the Euler flux
  // contribution below
  calcSAT_diff(params, roe_vars, roe_vars_dot,  dq, nrm, fluxL_dot, fluxR_dot);

//  calcSAT(params, nrm, dq, sat, u, v, H, use_efix)
  
  euler_fluxjac = params.euler_fluxjac;
  Tsol euler_fluxjac[numDofPerNode*numDofPerNode];

  Tsol v_vals[MAX_NUMDOFPERNODE];
//  nrm2 = params.nrm
  Tmsh nrm2[2];
  nrm2[1] = nx;   // why are we assigning to nrm2?
  nrm2[2] = ny;

//  convertFromNaturalToWorkingVars(params, q, v_vals)
  calcEulerFlux_diff(params, q, aux_vars, nrm2, euler_fluxjac);

  for (int i=0; i < numDofPerNode*numDofPerNode; ++i)
    fluxL_dot[i] += euler_fluxjac[i];

} // function RoeSolver_diff


template <typename Tsol, typename Tmsh, typename Tres>
void calcSAT_diff(Tsol roe_vars[3], Tsol roe_vars_dot[16], Tsol* dq, 
                  Tmsh* nrm, Tres sat_jacL*, Tres sat_jacR*, int use_efix=1)
{

// roe_vars = [u, v, H] at Roe average 
// roe_vars_dot contains all the non-zero derivatives of the roe_vars packed
// into a vector

  const int numDofPerNode = 4;
  // dq_dotL* = 1, dq_dotR* = -1, so don't pass them explicitly
  // SAT parameters
  double sat_Vn = 0.025;
  double sat_Vl = 0.025;
  double tau = 1.0;

  auto u = roe_vars[1];
  auto v = roe_vars[2];
  auto H = roe_vars[3];

  auto u_dotL1 = roe_vars_dot[1];
  auto u_dotR1 = roe_vars_dot[2];
  auto u_dotL2 = roe_vars_dot[3];
  auto u_dotR2 = roe_vars_dot[4];

  auto v_dotL1 = roe_vars_dot[5];
  auto v_dotR1 = roe_vars_dot[6];
  auto v_dotL3 = roe_vars_dot[7];
  auto v_dotR3 = roe_vars_dot[8];

  auto H_dotL1 = roe_vars_dot[9];
  auto H_dotR1 = roe_vars_dot[10];
  auto H_dotL2 = roe_vars_dot[11];
  auto H_dotR2 = roe_vars_dot[12];
  auto H_dotL3 = roe_vars_dot[13];
  auto H_dotR3 = roe_vars_dot[14];
  auto H_dotL4 = roe_vars_dot[15];
  auto H_dotR4 = roe_vars_dot[16];

//  @printit u v H

//  @printit u_dotL1 u_dotR1 u_dotL2 u_dotR2 v_dotL1 v_dotR1 v_dotL3 v_dotR3 H_dotL1 H_dotR1 H_dotL2 H_dotR2 H_dotL3 H_dotR3 H_dotL4 H_dotR4

  const double gami = 0.4;

  // Begin main execution
  auto nx = nrm[1];
  auto ny = nrm[2];

  auto dA = std::sqrt(nx*nx + ny*ny);

  auto Un = u*nx + v*ny // Normal Velocity;
  auto Un_dotL1 = nx*u_dotL1 + ny*v_dotL1;
  auto Un_dotR1 = nx*u_dotR1 + ny*v_dotR1;

  auto Un_dotL2 = nx*u_dotL2;
  auto Un_dotR2 = nx*u_dotR2;

  auto Un_dotL3 = ny*v_dotL3;
  auto Un_dotR3 = ny*v_dotR3;

  auto phi = 0.5*(u*u + v*v);

  auto phi_dotL1 = u*u_dotL1 + v*v_dotL1;
  auto phi_dotR1 = u*u_dotR1 + v*v_dotR1;

  auto phi_dotL2 = u*u_dotL2;
  auto phi_dotR2 = u*u_dotR2;
  
  auto phi_dotL3 = v*v_dotL3;
  auto phi_dotR3 = v*v_dotR3;


//  @printit Un_dotL1 phi_dotL1

  auto a = std::sqrt(gami*(H - phi)); // speed of sound
  auto t1 = gami/(2*a);
  auto a_dotL1 = t1*(H_dotL1 - phi_dotL1);
  auto a_dotR1 = t1*(H_dotR1 - phi_dotR1);

  auto a_dotL2 = t1*(H_dotL2 - phi_dotL2);
  auto a_dotR2 = t1*(H_dotR2 - phi_dotR2);

  auto a_dotL3 = t1*(H_dotL3 - phi_dotL3);
  auto a_dotR3 = t1*(H_dotR3 - phi_dotR3);

  auto a_dotL4 = t1*H_dotL4;
  auto a_dotR4 = t1*H_dotR4;

//  @printit a_dotL1

  auto lambda1 = Un + dA*a;
  auto lambda1_dotL1 = Un_dotL1 + dA*a_dotL1;
  auto lambda1_dotR1 = Un_dotR1 + dA*a_dotR1;

  auto lambda1_dotL2 = Un_dotL2 + dA*a_dotL2;
  auto lambda1_dotR2 = Un_dotR2 + dA*a_dotR2;

  auto lambda1_dotL3 = Un_dotL3 + dA*a_dotL3;
  auto lambda1_dotR3 = Un_dotR3 + dA*a_dotR3;

  auto lambda1_dotL4 = dA*a_dotL4;
  auto lambda1_dotR4 = dA*a_dotR4;

  
  auto lambda2 = Un - dA*a;
  auto lambda2_dotL1 = Un_dotL1 - dA*a_dotL1;
  auto lambda2_dotR1 = Un_dotR1 - dA*a_dotR1;

  auto lambda2_dotL2 = Un_dotL2 - dA*a_dotL2;
  auto lambda2_dotR2 = Un_dotR2 - dA*a_dotR2;

  auto lambda2_dotL3 = Un_dotL3 - dA*a_dotL3;
  auto lambda2_dotR3 = Un_dotR3 - dA*a_dotR3;

  auto lambda2_dotL4 = -dA*a_dotL4;
  auto lambda2_dotR4 = -dA*a_dotR4;


  auto lambda3 = Un;
  auto lambda3_dotL1 = Un_dotL1;
  auto lambda3_dotR1 = Un_dotR1;

  auto lambda3_dotL2 = Un_dotL2;
  auto lambda3_dotR2 = Un_dotR2;

  auto lambda3_dotL3 = Un_dotL3;
  auto lambda3_dotR3 = Un_dotR3;

//  @printit lambda1 lambda2 lambda3 lambda1_dotR1 lambda2_dotR1 lambda3_dotR1

  auto rhoA = absvalue(Un) + dA*a;
  //TODO: see if there is a better way to do this
  if (Un > 0)
    fac = 1;
  else
    fac = -1;

  auto rhoA_dotL1 = fac*Un_dotL1 + dA*a_dotL1;
  auto rhoA_dotR1 = fac*Un_dotR1 + dA*a_dotR1;

  auto rhoA_dotL2 = fac*Un_dotL2 + dA*a_dotL2;
  auto rhoA_dotR2 = fac*Un_dotR2 + dA*a_dotR2;

  auto rhoA_dotL3 = fac*Un_dotL3 + dA*a_dotL3;
  auto rhoA_dotR3 = fac*Un_dotR3 + dA*a_dotR3;

  auto rhoA_dotL4 = dA*a_dotL4;
  auto rhoA_dotR4 = dA*a_dotR4;

//  @printit rhoA rhoA_dotR1 sat_Vn

  // Compute Eigen Values of the Flux Jacobian
  // The eigen values calculated above cannot be used directly. Near stagnation
  // points lambda3 approaches zero while near sonic lines lambda1 and lambda2
  // approach zero. This has a possibility of creating numerical difficulties.
  // As a result, the eigen values are limited by the following expressions.


  //TODO: this does not handle use_efix = 0 case correctly

  // see lambda1 expression below
  if (absvalue(lambda1) > sat_Vn*rhoA)
  {
//    println("lambda1 is used")
    if (lambda1 > 0)
      fac = 1;
    else
      fac = -1;

//    @printit fac tau lambda1 lambda2 lambda3
//    println("fac = ", fac)

    t1 = tau*fac;
    //TODO: lambda1_dotL1 - lambgda1_dotL1 = 0, so simplify this
    auto lambda1_dotL1 = use_efix * 0.5 * (t1 * lambda1_dotL1 - lambda1_dotL1) + (1-use_efix)*lambda1_dotL1;
    auto lambda1_dotR1 = use_efix * 0.5 * (t1 * lambda1_dotR1 - lambda1_dotR1) + (1-use_efix)*lambda1_dotR1;

    auto lambda1_dotL2 = use_efix * 0.5 * (t1 * lambda1_dotL2 - lambda1_dotL2) + (1-use_efix)*lambda1_dotL2;
    auto lambda1_dotR2 = use_efix * 0.5 * (t1 * lambda1_dotR2 - lambda1_dotR2) + (1-use_efix)*lambda1_dotR2;

    auto lambda1_dotL3 = use_efix * 0.5 * (t1 * lambda1_dotL3 - lambda1_dotL3) + (1-use_efix)*lambda1_dotL3;
    auto lambda1_dotR3 = use_efix * 0.5 * (t1 * lambda1_dotR3 - lambda1_dotR3) + (1-use_efix)*lambda1_dotR3;

    auto lambda1_dotL4 = use_efix * 0.5 * (t1 * lambda1_dotL4 - lambda1_dotL4) + (1-use_efix)*lambda1_dotL4;
    auto lambda1_dotR4 = use_efix * 0.5 * (t1 * lambda1_dotR4 - lambda1_dotR4) + (1-use_efix)*lambda1_dotR4;


  } else
  {
//    println("not using lambda1")
    t1 = sat_Vn*tau;
    auto lambda1_dotL1 = use_efix * 0.5 * (t1 * rhoA_dotL1 - lambda1_dotL1) + (1-use_efix)*lambda1_dotL1;
    auto lambda1_dotR1 = use_efix * 0.5 * (t1 * rhoA_dotR1 - lambda1_dotR1) + (1-use_efix)*lambda1_dotR1;
 
    auto lambda1_dotL2 = use_efix * 0.5 * (t1 * rhoA_dotL2 - lambda1_dotL2) + (1-use_efix)*lambda1_dotL2;
    auto lambda1_dotR2 = use_efix * 0.5 * (t1 * rhoA_dotR2 - lambda1_dotR2) + (1-use_efix)*lambda1_dotR2;
    
    auto lambda1_dotL3 = use_efix * 0.5 * (t1 * rhoA_dotL3 - lambda1_dotL3) + (1-use_efix)*lambda1_dotL3;
    auto lambda1_dotR3 = use_efix * 0.5 * (t1 * rhoA_dotR3 - lambda1_dotR3) + (1-use_efix)*lambda1_dotR3;
 
    auto lambda1_dotL4 = use_efix * 0.5 * (t1 * rhoA_dotL4 - lambda1_dotL4) + (1-use_efix)*lambda1_dotL4;
    auto lambda1_dotR4 = use_efix * 0.5 * (t1 * rhoA_dotR4 - lambda1_dotR4) + (1-use_efix)*lambda1_dotR4;
 
  }

  auto lambda1 = use_efix*0.5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) - lambda1) + (1-use_efix)*lambda1;

  // see lambda2 expression below
  if (absvalue(lambda2) > sat_Vn*rhoA)
  {
    if (lambda2 > 0)
      fac = 1;
    else
      fac = -1;

    t1 = tau*fac;
    auto lambda2_dotL1 = use_efix * 0.5 * (t1 * lambda2_dotL1 - lambda2_dotL1) + (1-use_efix)*lambda2_dotL1;
    auto lambda2_dotR1 = use_efix * 0.5 * (t1 * lambda2_dotR1 - lambda2_dotR1) + (1-use_efix)*lambda2_dotR1;

    auto lambda2_dotL2 = use_efix * 0.5 * (t1 * lambda2_dotL2 - lambda2_dotL2) + (1-use_efix)*lambda2_dotL2;
    auto lambda2_dotR2 = use_efix * 0.5 * (t1 * lambda2_dotR2 - lambda2_dotR2) + (1-use_efix)*lambda2_dotR2;

    auto lambda2_dotL3 = use_efix * 0.5 * (t1 * lambda2_dotL3 - lambda2_dotL3) + (1-use_efix)*lambda2_dotL3;
    auto lambda2_dotR3 = use_efix * 0.5 * (t1 * lambda2_dotR3 - lambda2_dotR3) + (1-use_efix)*lambda2_dotR3;

    auto lambda2_dotL4 = use_efix * 0.5 * (t1 * lambda2_dotL4 - lambda2_dotL4) + (1-use_efix)*lambda2_dotL4;
    auto lambda2_dotR4 = use_efix * 0.5 * (t1 * lambda2_dotR4 - lambda2_dotR4) + (1-use_efix)*lambda2_dotR4;

  } else
  {
    t1 = sat_Vn*tau;
    auto lambda2_dotL1 = use_efix * 0.5 * (t1 * rhoA_dotL1 - lambda2_dotL1) + (1-use_efix)*lambda2_dotL1;
    auto lambda2_dotR1 = use_efix * 0.5 * (t1 * rhoA_dotR1 - lambda2_dotR1) + (1-use_efix)*lambda2_dotR1;
 
    auto lambda2_dotL2 = use_efix * 0.5 * (t1 * rhoA_dotL2 - lambda2_dotL2) + (1-use_efix)*lambda2_dotL2;
    auto lambda2_dotR2 = use_efix * 0.5 * (t1 * rhoA_dotR2 - lambda2_dotR2) + (1-use_efix)*lambda2_dotR2;
    
    auto lambda2_dotL3 = use_efix * 0.5 * (t1 * rhoA_dotL3 - lambda2_dotL3) + (1-use_efix)*lambda2_dotL3;
    auto lambda2_dotR3 = use_efix * 0.5 * (t1 * rhoA_dotR3 - lambda2_dotR3) + (1-use_efix)*lambda2_dotR3;
 
    auto lambda2_dotL4 = use_efix * 0.5 * (t1 * rhoA_dotL4 - lambda2_dotL4) + (1-use_efix)*lambda2_dotL4;
    auto lambda2_dotR4 = use_efix * 0.5 * (t1 * rhoA_dotR4 - lambda2_dotR4) + (1-use_efix)*lambda2_dotR4;
 
  }

  auto lambda2 = use_efix*0.5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) - lambda2) + (1-use_efix)*lambda2;


  // see lambda3 expression below
  if (absvalue(lambda3) > sat_Vn*rhoA)
  {
    if (lambda3 > 0)
      fac = 1;
    else
      fac = -1;

    t1 = tau*fac;
    auto lambda3_dotL1 = use_efix * 0.5 * (t1 * lambda3_dotL1 - lambda3_dotL1) + (1-use_efix)*lambda3_dotL1;
    auto lambda3_dotR1 = use_efix * 0.5 * (t1 * lambda3_dotR1 - lambda3_dotR1) + (1-use_efix)*lambda3_dotR1;

    auto lambda3_dotL2 = use_efix * 0.5 * (t1 * lambda3_dotL2 - lambda3_dotL2) + (1-use_efix)*lambda3_dotL2;
    auto lambda3_dotR2 = use_efix * 0.5 * (t1 * lambda3_dotR2 - lambda3_dotR2) + (1-use_efix)*lambda3_dotR2;

    auto lambda3_dotL3 = use_efix * 0.5 * (t1 * lambda3_dotL3 - lambda3_dotL3) + (1-use_efix)*lambda3_dotL3;
    auto lambda3_dotR3 = use_efix * 0.5 * (t1 * lambda3_dotR3 - lambda3_dotR3) + (1-use_efix)*lambda3_dotR3;

    auto lambda3_dotL4 = 0.0;
    auto lambda3_dotR4 = 0.0;


  } else
  {
    t1 = sat_Vn*tau;
    auto lambda3_dotL1 = use_efix * 0.5 * (t1 * rhoA_dotL1 - lambda3_dotL1) + (1-use_efix)*lambda3_dotL1;
    auto lambda3_dotR1 = use_efix * 0.5 * (t1 * rhoA_dotR1 - lambda3_dotR1) + (1-use_efix)*lambda3_dotR1;
 
    auto lambda3_dotL2 = use_efix * 0.5 * (t1 * rhoA_dotL2 - lambda3_dotL2) + (1-use_efix)*lambda3_dotL2;
    auto lambda3_dotR2 = use_efix * 0.5 * (t1 * rhoA_dotR2 - lambda3_dotR2) + (1-use_efix)*lambda3_dotR2;
    
    auto lambda3_dotL3 = use_efix * 0.5 * (t1 * rhoA_dotL3 - lambda3_dotL3) + (1-use_efix)*lambda3_dotL3;
    auto lambda3_dotR3 = use_efix * 0.5 * (t1 * rhoA_dotR3 - lambda3_dotR3) + (1-use_efix)*lambda3_dotR3;
 

    auto lambda3_dotL4 = use_efix * 0.5 * t1 * rhoA_dotL4 + (1-use_efix)*0.0;
    auto lambda3_dotR4 = use_efix * 0.5 * t1 * rhoA_dotR4 + (1-use_efix)*0.0;
 
  }

  auto lambda3 = use_efix*0.5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) - lambda3) + (1-use_efix)*lambda3;

//  println("after entropy fix")
//  @printit lambda1 lambda2 lambda3 lambda1_dotR1 lambda2_dotR1 lambda3_dotR1

                    
  auto dq1 = dq[1];
  auto dq2 = dq[2];
  auto dq3 = dq[3];
  auto dq4 = dq[4];

  //TODO: see if the sat values are needed
  // sat[1] = lambda3*dq1
  sat_jacL[1, 1] = lambda3 + dq1*lambda3_dotL1;
  sat_jacL[1, 2] =         + dq1*lambda3_dotL2;
  sat_jacL[1, 3] =         + dq1*lambda3_dotL3;
  sat_jacL[1, 4] =         + dq1*lambda3_dotL4;
  //TODO: zero out unused entries?

  sat_jacR[1, 1] = -lambda3 + dq1*lambda3_dotR1;
  sat_jacR[1, 2] =          + dq1*lambda3_dotR2;
  sat_jacR[1, 3] =          + dq1*lambda3_dotR3;
  sat_jacR[1, 4] =          + dq1*lambda3_dotR4;

  // sat[2] = lambda3*dq2
  sat_jacL[2, 1] =         + dq2*lambda3_dotL1;
  sat_jacL[2, 2] = lambda3 + dq2*lambda3_dotL2;
  sat_jacL[2, 3] =         + dq2*lambda3_dotL3;
  sat_jacL[2, 4] =         + dq2*lambda3_dotL4;

  sat_jacR[2, 1] =          + dq2*lambda3_dotR1;
  sat_jacR[2, 2] = -lambda3 + dq2*lambda3_dotR2;
  sat_jacR[2, 3] =          + dq2*lambda3_dotR3;
  sat_jacR[2, 4] =          + dq2*lambda3_dotR4;

  // sat[3] = lambda3*dq3
  sat_jacL[3, 1] =         + dq3*lambda3_dotL1;
  sat_jacL[3, 2] =         + dq3*lambda3_dotL2;
  sat_jacL[3, 3] = lambda3 + dq3*lambda3_dotL3;
  sat_jacL[3, 4] =         + dq3*lambda3_dotL4;

  sat_jacR[3, 1] =          + dq3*lambda3_dotR1;
  sat_jacR[3, 2] =          + dq3*lambda3_dotR2;
  sat_jacR[3, 3] = -lambda3 + dq3*lambda3_dotR3;
  sat_jacR[3, 4] =          + dq3*lambda3_dotR4;

  // sat[4] = lambda3*dq4;
  sat_jacL[4, 1] =           dq4*lambda3_dotL1;
  sat_jacL[4, 2] =           dq4*lambda3_dotL2;
  sat_jacL[4, 3] =           dq4*lambda3_dotL3;
  sat_jacL[4, 4] = lambda3 + dq4*lambda3_dotL4;

  sat_jacR[4, 1] =            dq4*lambda3_dotR1;
  sat_jacR[4, 2] =            dq4*lambda3_dotR2;
  sat_jacR[4, 3] =            dq4*lambda3_dotR3;
  sat_jacR[4, 4] = -lambda3 + dq4*lambda3_dotR4;

//  @printit sat_jacR[2, 1]

//  E1dq = params.res_vals1
//  E2dq = params.res_vals2

  Tres E1dq[MAX_DOFPERNODE];
  Tres E2dq[MAX_DOFPERNODE];

  //-- get E1*dq
  E1dq[1] = phi*dq1 - u*dq2 - v*dq3 + dq4;
  auto E1dq1_dotL1 = phi + dq1*phi_dotL1     - dq2*u_dotL1     - dq3*v_dotL1;
  auto E1dq1_dotL2 =       dq1*phi_dotL2 - u - dq2*u_dotL2;
  auto E1dq1_dotL3 =       dq1*phi_dotL3                   - v - dq3*v_dotL3;
  auto E1dq1_dotL4 = 1;

  auto E1dq1_dotR1 = -phi + dq1*phi_dotR1     - dq2*u_dotR1     - dq3*v_dotR1;
  auto E1dq1_dotR2 =        dq1*phi_dotR2 + u - dq2*u_dotR2;
  auto E1dq1_dotR3 =        dq1*phi_dotR3                   + v - dq3*v_dotR3;
  auto E1dq1_dotR4 = -1;


  E1dq[2] = E1dq[1]*u;
  auto E1dq2_dotL1 = u*E1dq1_dotL1 + E1dq[1]*u_dotL1;
  auto E1dq2_dotL2 = u*E1dq1_dotL2 + E1dq[1]*u_dotL2;
  auto E1dq2_dotL3 = u*E1dq1_dotL3;
  auto E1dq2_dotL4 = u*E1dq1_dotL4;

  auto E1dq2_dotR1 = u*E1dq1_dotR1 + E1dq[1]*u_dotR1;
  auto E1dq2_dotR2 = u*E1dq1_dotR2 + E1dq[1]*u_dotR2;
  auto E1dq2_dotR3 = u*E1dq1_dotR3;
  auto E1dq2_dotR4 = u*E1dq1_dotR4;

  E1dq[3] = E1dq[1]*v;
  auto E1dq3_dotL1 = v*E1dq1_dotL1 + E1dq[1]*v_dotL1;
  auto E1dq3_dotL2 = v*E1dq1_dotL2;
  auto E1dq3_dotL3 = v*E1dq1_dotL3 + E1dq[1]*v_dotL3;
  auto E1dq3_dotL4 = v*E1dq1_dotL4;

  auto E1dq3_dotR1 = v*E1dq1_dotR1 + E1dq[1]*v_dotR1;
  auto E1dq3_dotR2 = v*E1dq1_dotR2;
  auto E1dq3_dotR3 = v*E1dq1_dotR3 + E1dq[1]*v_dotR3;
  auto E1dq3_dotR4 = v*E1dq1_dotR4;

  E1dq[4] = E1dq[1]*H;
  auto E1dq4_dotL1 = H*E1dq1_dotL1 + E1dq[1]*H_dotL1;
  auto E1dq4_dotL2 = H*E1dq1_dotL2 + E1dq[1]*H_dotL2;
  auto E1dq4_dotL3 = H*E1dq1_dotL3 + E1dq[1]*H_dotL3;
  auto E1dq4_dotL4 = H*E1dq1_dotL4 + E1dq[1]*H_dotL4;

  auto E1dq4_dotR1 = H*E1dq1_dotR1 + E1dq[1]*H_dotR1;
  auto E1dq4_dotR2 = H*E1dq1_dotR2 + E1dq[1]*H_dotR2;
  auto E1dq4_dotR3 = H*E1dq1_dotR3 + E1dq[1]*H_dotR3;
  auto E1dq4_dotR4 = H*E1dq1_dotR4 + E1dq[1]*H_dotR4;



  //-- get E2*dq
  E2dq[1] = 0.0;
  E2dq[2] = -Un*dq1 + nx*dq2 + ny*dq3;
  auto E2dq2_dotL1 = -Un + -Un_dotL1*dq1;
  auto E2dq2_dotL2 =     + -Un_dotL2*dq1 + nx;
  auto E2dq2_dotL3 =     + -Un_dotL3*dq1      + ny;

  auto E2dq2_dotR1 = Un + -Un_dotR1*dq1;
  auto E2dq2_dotR2 =    + -Un_dotR2*dq1 - nx;
  auto E2dq2_dotR3 =    + -Un_dotR3*dq1      - ny;


  E2dq[3] = E2dq[2]*ny;
  auto E2dq3_dotL1 = ny*E2dq2_dotL1;
  auto E2dq3_dotL2 = ny*E2dq2_dotL2;
  auto E2dq3_dotL3 = ny*E2dq2_dotL3;

  auto E2dq3_dotR1 = ny*E2dq2_dotR1;
  auto E2dq3_dotR2 = ny*E2dq2_dotR2;
  auto E2dq3_dotR3 = ny*E2dq2_dotR3;

  E2dq[4] = E2dq[2]*Un;
  auto E2dq4_dotL1 = Un*E2dq2_dotL1 + E2dq[2]*Un_dotL1;
  auto E2dq4_dotL2 = Un*E2dq2_dotL2 + E2dq[2]*Un_dotL2;
  auto E2dq4_dotL3 = Un*E2dq2_dotL3 + E2dq[2]*Un_dotL3;

  auto E2dq4_dotR1 = Un*E2dq2_dotR1 + E2dq[2]*Un_dotR1;
  auto E2dq4_dotR2 = Un*E2dq2_dotR2 + E2dq[2]*Un_dotR2;
  auto E2dq4_dotR3 = Un*E2dq2_dotR3 + E2dq[2]*Un_dotR3;

  E2dq[2] = E2dq[2]*nx;
  auto E2dq2_dotL1 = nx*E2dq2_dotL1;
  auto E2dq2_dotL2 = nx*E2dq2_dotL2;
  auto E2dq2_dotL3 = nx*E2dq2_dotL3;

  auto E2dq2_dotR1 = nx*E2dq2_dotR1;
  auto E2dq2_dotR2 = nx*E2dq2_dotR2;
  auto E2dq2_dotR3 = nx*E2dq2_dotR3;

//  @printit E1dq1_dotR1 E1dq2_dotR1 E1dq3_dotR1 E1dq4_dotR1 E2dq2_dotR1

  //-- add to sat
  auto tmp1 = 0.5*(lambda1 + lambda2) - lambda3;
  auto tmp1_dotL1 = 0.5*(lambda1_dotL1 + lambda2_dotL1) - lambda3_dotL1;
  auto tmp1_dotL2 = 0.5*(lambda1_dotL2 + lambda2_dotL2) - lambda3_dotL2;
  auto tmp1_dotL3 = 0.5*(lambda1_dotL3 + lambda2_dotL3) - lambda3_dotL3;
  auto tmp1_dotL4 = 0.5*(lambda1_dotL4 + lambda2_dotL4) - lambda3_dotL4;

  auto tmp1_dotR1 = 0.5*(lambda1_dotR1 + lambda2_dotR1) - lambda3_dotR1;
  auto tmp1_dotR2 = 0.5*(lambda1_dotR2 + lambda2_dotR2) - lambda3_dotR2;
  auto tmp1_dotR3 = 0.5*(lambda1_dotR3 + lambda2_dotR3) - lambda3_dotR3;
  auto tmp1_dotR4 = 0.5*(lambda1_dotR4 + lambda2_dotR4) - lambda3_dotR4;


  auto tmp2 = gami/(a*a);
  t1 = -2*tmp2/a;
//  t1 = -2*gami/(a*a*a) // = -2*tmp2/a;
  auto tmp2_dotL1 = t1*a_dotL1;
  auto tmp2_dotL2 = t1*a_dotL2;
  auto tmp2_dotL3 = t1*a_dotL3;
  auto tmp2_dotL4 = t1*a_dotL4;

  auto tmp2_dotR1 = t1*a_dotR1;
  auto tmp2_dotR2 = t1*a_dotR2;
  auto tmp2_dotR3 = t1*a_dotR3;
  auto tmp2_dotR4 = t1*a_dotR4;

  auto tmp3 = 1.0/(dA*dA);
  //for i=1:length(sat)
  //  sat[i] = sat[i] + tmp1*(tmp2*E1dq[i] + tmp3*E2dq[i])
  //end

//  @printit tmp1 tmp2 E1dq1_dotR1 E1dq[1] tmp1_dotR1 tmp2_dotR1 tmp3 E2dq[1]

  //TODO: align + signs
  sat_jacL[1,1] += tmp1*tmp2*E1dq1_dotL1   + tmp1*E1dq[1]*tmp2_dotL1 + 
                   tmp2*E1dq[1]*tmp1_dotL1 + 
                                           + tmp3*E2dq[1]*tmp1_dotL1;

  sat_jacL[1,2] += tmp1*tmp2*E1dq1_dotL2   + tmp1*E1dq[1]*tmp2_dotL2 + 
                   tmp2*E1dq[1]*tmp1_dotL2 + 
                                           + tmp3*E2dq[1]*tmp1_dotL2;

  sat_jacL[1,3] += tmp1*tmp2*E1dq1_dotL3   + tmp1*E1dq[1]*tmp2_dotL3 + 
                   tmp2*E1dq[1]*tmp1_dotL3 + 
                                           + tmp3*E2dq[1]*tmp1_dotL3;

  sat_jacL[1,4] += tmp1*tmp2*E1dq1_dotL4   + tmp1*E1dq[1]*tmp2_dotL4 + 
                   tmp2*E1dq[1]*tmp1_dotL4 + 
                                           + tmp3*E2dq[1]*tmp1_dotL4;

  sat_jacR[1,1] += tmp1*tmp2*E1dq1_dotR1   + tmp1*E1dq[1]*tmp2_dotR1 + 
                   tmp2*E1dq[1]*tmp1_dotR1 + 
                                           + tmp3*E2dq[1]*tmp1_dotR1;

  sat_jacR[1,2] += tmp1*tmp2*E1dq1_dotR2   + tmp1*E1dq[1]*tmp2_dotR2 + 
                   tmp2*E1dq[1]*tmp1_dotR2 + 
                                           + tmp3*E2dq[1]*tmp1_dotR2;

  sat_jacR[1,3] += tmp1*tmp2*E1dq1_dotR3   + tmp1*E1dq[1]*tmp2_dotR3 + 
                   tmp2*E1dq[1]*tmp1_dotR3 + 
                                           + tmp3*E2dq[1]*tmp1_dotR3;

  sat_jacR[1,4] += tmp1*tmp2*E1dq1_dotR4   + tmp1*E1dq[1]*tmp2_dotR4 + 
                   tmp2*E1dq[1]*tmp1_dotR4 + 
                                           + tmp3*E2dq[1]*tmp1_dotR4;


  sat_jacL[2,1] += tmp1*tmp2*E1dq2_dotL1   + tmp1*E1dq[2]*tmp2_dotL1 + 
                   tmp2*E1dq[2]*tmp1_dotL1 + 
                   tmp1*tmp3*E2dq2_dotL1   + tmp3*E2dq[2]*tmp1_dotL1;

  sat_jacL[2,2] += tmp1*tmp2*E1dq2_dotL2   + tmp1*E1dq[2]*tmp2_dotL2 + 
                   tmp2*E1dq[2]*tmp1_dotL2 + 
                   tmp1*tmp3*E2dq2_dotL2   + tmp3*E2dq[2]*tmp1_dotL2;

  sat_jacL[2,3] += tmp1*tmp2*E1dq2_dotL3   + tmp1*E1dq[2]*tmp2_dotL3 + 
                   tmp2*E1dq[2]*tmp1_dotL3 + 
                   tmp1*tmp3*E2dq2_dotL3   + tmp3*E2dq[2]*tmp1_dotL3;

  sat_jacL[2,4] += tmp1*tmp2*E1dq2_dotL4   + tmp1*E1dq[2]*tmp2_dotL4 + 
                   tmp2*E1dq[2]*tmp1_dotL4 + 
                                           + tmp3*E2dq[2]*tmp1_dotL4;

  sat_jacR[2,1] += tmp1*tmp2*E1dq2_dotR1   + tmp1*E1dq[2]*tmp2_dotR1 + 
                   tmp2*E1dq[2]*tmp1_dotR1 + 
                   tmp1*tmp3*E2dq2_dotR1   + tmp3*E2dq[2]*tmp1_dotR1;

  sat_jacR[2,2] += tmp1*tmp2*E1dq2_dotR2   + tmp1*E1dq[2]*tmp2_dotR2 + 
                   tmp2*E1dq[2]*tmp1_dotR2 + 
                   tmp1*tmp3*E2dq2_dotR2   + tmp3*E2dq[2]*tmp1_dotR2;

  sat_jacR[2,3] += tmp1*tmp2*E1dq2_dotR3   + tmp1*E1dq[2]*tmp2_dotR3 + 
                   tmp2*E1dq[2]*tmp1_dotR3 + 
                   tmp1*tmp3*E2dq2_dotR3   + tmp3*E2dq[2]*tmp1_dotR3;

  sat_jacR[2,4] += tmp1*tmp2*E1dq2_dotR4   + tmp1*E1dq[2]*tmp2_dotR4 + 
                   tmp2*E1dq[2]*tmp1_dotR4 + 
                                           + tmp3*E2dq[2]*tmp1_dotR4;


  sat_jacL[3,1] += tmp1*tmp2*E1dq3_dotL1   + tmp1*E1dq[3]*tmp2_dotL1 + 
                   tmp2*E1dq[3]*tmp1_dotL1 + 
                   tmp1*tmp3*E2dq3_dotL1   + tmp3*E2dq[3]*tmp1_dotL1;

  sat_jacL[3,2] += tmp1*tmp2*E1dq3_dotL2   + tmp1*E1dq[3]*tmp2_dotL2 + 
                   tmp2*E1dq[3]*tmp1_dotL2 + 
                   tmp1*tmp3*E2dq3_dotL2   + tmp3*E2dq[3]*tmp1_dotL2;

  sat_jacL[3,3] += tmp1*tmp2*E1dq3_dotL3   + tmp1*E1dq[3]*tmp2_dotL3 + 
                   tmp2*E1dq[3]*tmp1_dotL3 + 
                   tmp1*tmp3*E2dq3_dotL3   + tmp3*E2dq[3]*tmp1_dotL3;

  sat_jacL[3,4] += tmp1*tmp2*E1dq3_dotL4   + tmp1*E1dq[3]*tmp2_dotL4 + 
                   tmp2*E1dq[3]*tmp1_dotL4 + 
                                           + tmp3*E2dq[3]*tmp1_dotL4;

  sat_jacR[3,1] += tmp1*tmp2*E1dq3_dotR1   + tmp1*E1dq[3]*tmp2_dotR1 + 
                   tmp2*E1dq[3]*tmp1_dotR1 + 
                   tmp1*tmp3*E2dq3_dotR1   + tmp3*E2dq[3]*tmp1_dotR1;

  sat_jacR[3,2] += tmp1*tmp2*E1dq3_dotR2   + tmp1*E1dq[3]*tmp2_dotR2 + 
                   tmp2*E1dq[3]*tmp1_dotR2 + 
                   tmp1*tmp3*E2dq3_dotR2   + tmp3*E2dq[3]*tmp1_dotR2;

  sat_jacR[3,3] += tmp1*tmp2*E1dq3_dotR3   + tmp1*E1dq[3]*tmp2_dotR3 + 
                   tmp2*E1dq[3]*tmp1_dotR3 + 
                   tmp1*tmp3*E2dq3_dotR3   + tmp3*E2dq[3]*tmp1_dotR3;

  sat_jacR[3,4] += tmp1*tmp2*E1dq3_dotR4   + tmp1*E1dq[3]*tmp2_dotR4 + 
                   tmp2*E1dq[3]*tmp1_dotR4 + 
                                           + tmp3*E2dq[3]*tmp1_dotR4;


  sat_jacL[4,1] += tmp1*tmp2*E1dq4_dotL1   + tmp1*E1dq[4]*tmp2_dotL1 + 
                   tmp2*E1dq[4]*tmp1_dotL1 + 
                   tmp1*tmp3*E2dq4_dotL1   + tmp3*E2dq[4]*tmp1_dotL1;

  sat_jacL[4,2] += tmp1*tmp2*E1dq4_dotL2   + tmp1*E1dq[4]*tmp2_dotL2 + 
                   tmp2*E1dq[4]*tmp1_dotL2 + 
                   tmp1*tmp3*E2dq4_dotL2   + tmp3*E2dq[4]*tmp1_dotL2;

  sat_jacL[4,3] += tmp1*tmp2*E1dq4_dotL3   + tmp1*E1dq[4]*tmp2_dotL3 + 
                   tmp2*E1dq[4]*tmp1_dotL3 + 
                   tmp1*tmp3*E2dq4_dotL3   + tmp3*E2dq[4]*tmp1_dotL3;

  sat_jacL[4,4] += tmp1*tmp2*E1dq4_dotL4   + tmp1*E1dq[4]*tmp2_dotL4 + 
                   tmp2*E1dq[4]*tmp1_dotL4 + 
                                           + tmp3*E2dq[4]*tmp1_dotL4;

  sat_jacR[4,1] += tmp1*tmp2*E1dq4_dotR1   + tmp1*E1dq[4]*tmp2_dotR1 + 
                   tmp2*E1dq[4]*tmp1_dotR1 + 
                   tmp1*tmp3*E2dq4_dotR1   + tmp3*E2dq[4]*tmp1_dotR1;

  sat_jacR[4,2] += tmp1*tmp2*E1dq4_dotR2   + tmp1*E1dq[4]*tmp2_dotR2 + 
                   tmp2*E1dq[4]*tmp1_dotR2 + 
                   tmp1*tmp3*E2dq4_dotR2   + tmp3*E2dq[4]*tmp1_dotR2;

  sat_jacR[4,3] += tmp1*tmp2*E1dq4_dotR3   + tmp1*E1dq[4]*tmp2_dotR3 + 
                   tmp2*E1dq[4]*tmp1_dotR3 + 
                   tmp1*tmp3*E2dq4_dotR3   + tmp3*E2dq[4]*tmp1_dotR3;

  sat_jacR[4,4] += tmp1*tmp2*E1dq4_dotR4   + tmp1*E1dq[4]*tmp2_dotR4 + 
                   tmp2*E1dq[4]*tmp1_dotR4 + 
                                           + tmp3*E2dq[4]*tmp1_dotR4;

//  println("after first summation")
//  @printit sat_jacR[2, 1]

  //-- get E3*dq
  E1dq[1] = -Un*dq1 + nx*dq2 + ny*dq3;
  E1dq1_dotL1 = -Un - dq1*Un_dotL1;
  E1dq1_dotL2 =     - dq1*Un_dotL2 + nx;
  E1dq1_dotL3 =     - dq1*Un_dotL3      + ny;

  E1dq1_dotR1 = Un - dq1*Un_dotR1;
  E1dq1_dotR2 =    - dq1*Un_dotR2 - nx;
  E1dq1_dotR3 =    - dq1*Un_dotR3      - ny;


  E1dq[2] = E1dq[1]*u;
  E1dq2_dotL1 = u*E1dq1_dotL1 + E1dq[1]*u_dotL1;
  E1dq2_dotL2 = u*E1dq1_dotL2 + E1dq[1]*u_dotL2;
  E1dq2_dotL3 = u*E1dq1_dotL3;

  E1dq2_dotR1 = u*E1dq1_dotR1 + E1dq[1]*u_dotR1;
  E1dq2_dotR2 = u*E1dq1_dotR2 + E1dq[1]*u_dotR2;
  E1dq2_dotR3 = u*E1dq1_dotR3;


  E1dq[3] = E1dq[1]*v;
  E1dq3_dotL1 = v*E1dq1_dotL1 + E1dq[1]*v_dotL1;
  E1dq3_dotL2 = v*E1dq1_dotL2;
  E1dq3_dotL3 = v*E1dq1_dotL3 + E1dq[1]*v_dotL3;

  E1dq3_dotR1 = v*E1dq1_dotR1 + E1dq[1]*v_dotR1;
  E1dq3_dotR2 = v*E1dq1_dotR2;
  E1dq3_dotR3 = v*E1dq1_dotR3 + E1dq[1]*v_dotR3;


  E1dq[4] = E1dq[1]*H
  E1dq4_dotL1 = H*E1dq1_dotL1 + E1dq[1]*H_dotL1
  E1dq4_dotL2 = H*E1dq1_dotL2 + E1dq[1]*H_dotL2
  E1dq4_dotL3 = H*E1dq1_dotL3 + E1dq[1]*H_dotL3
  E1dq4_dotL4 =               + E1dq[1]*H_dotL4

  E1dq4_dotR1 = H*E1dq1_dotR1 + E1dq[1]*H_dotR1
  E1dq4_dotR2 = H*E1dq1_dotR2 + E1dq[1]*H_dotR2
  E1dq4_dotR3 = H*E1dq1_dotR3 + E1dq[1]*H_dotR3
  E1dq4_dotR4 =               + E1dq[1]*H_dotR4


  //-- get E4*dq
  t1 = phi*dq1 - u*dq2 - v*dq3 + dq4;
  t1_dotL1 = phi + dq1*phi_dotL1     - dq2*u_dotL1     - dq3*v_dotL1;
  t1_dotL2 =       dq1*phi_dotL2 - u - dq2*u_dotL2;
  t1_dotL3 =       dq1*phi_dotL3                   - v - dq3*v_dotL3;
  t1_dotL4 = 1;

  t1_dotR1 = -phi + dq1*phi_dotR1     - dq2*u_dotR1     - dq3*v_dotR1;
  t1_dotR2 =        dq1*phi_dotR2 + u - dq2*u_dotR2;
  t1_dotR3 =        dq1*phi_dotR3                   + v - dq3*v_dotR3;
  t1_dotR4 = -1;


  E2dq[1] = 0.0
  E2dq[2] = t1*nx
  E2dq2_dotL1 = nx*t1_dotL1;
  E2dq2_dotL2 = nx*t1_dotL2;
  E2dq2_dotL3 = nx*t1_dotL3;
  E2dq2_dotL4 = nx*t1_dotL4;

  E2dq2_dotR1 = nx*t1_dotR1;
  E2dq2_dotR2 = nx*t1_dotR2;
  E2dq2_dotR3 = nx*t1_dotR3;
  E2dq2_dotR4 = nx*t1_dotR4;


  E2dq[3] = t1*ny;
  E2dq3_dotL1 = ny*t1_dotL1;
  E2dq3_dotL2 = ny*t1_dotL2;
  E2dq3_dotL3 = ny*t1_dotL3;
  E2dq3_dotL4 = ny*t1_dotL4;

  E2dq3_dotR1 = ny*t1_dotR1;
  E2dq3_dotR2 = ny*t1_dotR2;
  E2dq3_dotR3 = ny*t1_dotR3;
  E2dq3_dotR4 = ny*t1_dotR4;


  E2dq[4] = t1*Un;
  E2dq4_dotL1 = Un*t1_dotL1 + t1*Un_dotL1;
  E2dq4_dotL2 = Un*t1_dotL2 + t1*Un_dotL2;
  E2dq4_dotL3 = Un*t1_dotL3 + t1*Un_dotL3;
  E2dq4_dotL4 = Un*t1_dotL4;

  E2dq4_dotR1 = Un*t1_dotR1 + t1*Un_dotR1;
  E2dq4_dotR2 = Un*t1_dotR2 + t1*Un_dotR2;
  E2dq4_dotR3 = Un*t1_dotR3 + t1*Un_dotR3;
  E2dq4_dotR4 = Un*t1_dotR4;

//  @printit E1dq1_dotR1 E1dq2_dotR1 E1dq3_dotR1 E1dq4_dotR1 E2dq2_dotR1 E2dq3_dotR1 E2dq4_dotR1

  //-- add to sat
  t1 = 1/(dA*a);
  t2 = 1/a;
  tmp1 = 0.5*(lambda1 - lambda2)*t1;

  tmp1_dotL1 = 0.5 * t1 * ( (lambda1_dotL1 - lambda2_dotL1) - (lambda1 - lambda2) * t2 * a_dotL1);
  tmp1_dotL2 = 0.5 * t1 * ( (lambda1_dotL2 - lambda2_dotL2) - (lambda1 - lambda2) * t2 * a_dotL2);
  tmp1_dotL3 = 0.5 * t1 * ( (lambda1_dotL3 - lambda2_dotL3) - (lambda1 - lambda2) * t2 * a_dotL3);
  tmp1_dotL4 = 0.5 * t1 * ( (lambda1_dotL4 - lambda2_dotL4) - (lambda1 - lambda2) * t2 * a_dotL4);

  tmp1_dotR1 = 0.5 * t1 * ( (lambda1_dotR1 - lambda2_dotR1) - (lambda1 - lambda2) * t2 * a_dotR1);
  tmp1_dotR2 = 0.5 * t1 * ( (lambda1_dotR2 - lambda2_dotR2) - (lambda1 - lambda2) * t2 * a_dotR2);
  tmp1_dotR3 = 0.5 * t1 * ( (lambda1_dotR3 - lambda2_dotR3) - (lambda1 - lambda2) * t2 * a_dotR3);
  tmp1_dotR4 = 0.5 * t1 * ( (lambda1_dotR4 - lambda2_dotR4) - (lambda1 - lambda2) * t2 * a_dotR4);


  //for i=1:length(sat)
  //  sat[i] = sat[i] + tmp1*(E1dq[i] + gami*E2dq[i])
  //end

  t1 = E1dq[1] + gami*E2dq[1]
  sat_jacL[1, 1] += tmp1*E1dq1_dotL1 + t1*tmp1_dotL1;
  sat_jacL[1, 2] += tmp1*E1dq1_dotL2 + t1*tmp1_dotL2;
  sat_jacL[1, 3] += tmp1*E1dq1_dotL3 + t1*tmp1_dotL3;
  sat_jacL[1, 4] +=                    t1*tmp1_dotL4;

  sat_jacR[1, 1] += tmp1*E1dq1_dotR1 + t1*tmp1_dotR1;
  sat_jacR[1, 2] += tmp1*E1dq1_dotR2 + t1*tmp1_dotR2;
  sat_jacR[1, 3] += tmp1*E1dq1_dotR3 + t1*tmp1_dotR3;
  sat_jacR[1, 4] +=                    t1*tmp1_dotR4;

  t1 = E1dq[2] + gami*E2dq[2];
  sat_jacL[2, 1] += tmp1*(E1dq2_dotL1 + gami*E2dq2_dotL1) + t1*tmp1_dotL1;
  sat_jacL[2, 2] += tmp1*(E1dq2_dotL2 + gami*E2dq2_dotL2) + t1*tmp1_dotL2;
  sat_jacL[2, 3] += tmp1*(E1dq2_dotL3 + gami*E2dq2_dotL3) + t1*tmp1_dotL3;
  sat_jacL[2, 4] += tmp1*(            + gami*E2dq2_dotL4) + t1*tmp1_dotL4;

  sat_jacR[2, 1] += tmp1*(E1dq2_dotR1 + gami*E2dq2_dotR1) + t1*tmp1_dotR1;
  sat_jacR[2, 2] += tmp1*(E1dq2_dotR2 + gami*E2dq2_dotR2) + t1*tmp1_dotR2;
  sat_jacR[2, 3] += tmp1*(E1dq2_dotR3 + gami*E2dq2_dotR3) + t1*tmp1_dotR3;
  sat_jacR[2, 4] += tmp1*(            + gami*E2dq2_dotR4) + t1*tmp1_dotR4;

  t1 = E1dq[3] + gami*E2dq[3];
  sat_jacL[3, 1] += tmp1*(E1dq3_dotL1 + gami*E2dq3_dotL1) + t1*tmp1_dotL1;
  sat_jacL[3, 2] += tmp1*(E1dq3_dotL2 + gami*E2dq3_dotL2) + t1*tmp1_dotL2;
  sat_jacL[3, 3] += tmp1*(E1dq3_dotL3 + gami*E2dq3_dotL3) + t1*tmp1_dotL3;
  sat_jacL[3, 4] += tmp1*(            + gami*E2dq3_dotL4) + t1*tmp1_dotL4;

  sat_jacR[3, 1] += tmp1*(E1dq3_dotR1 + gami*E2dq3_dotR1) + t1*tmp1_dotR1;
  sat_jacR[3, 2] += tmp1*(E1dq3_dotR2 + gami*E2dq3_dotR2) + t1*tmp1_dotR2;
  sat_jacR[3, 3] += tmp1*(E1dq3_dotR3 + gami*E2dq3_dotR3) + t1*tmp1_dotR3;
  sat_jacR[3, 4] += tmp1*(            + gami*E2dq3_dotR4) + t1*tmp1_dotR4;

  t1 = E1dq[4] + gami*E2dq[4];
  sat_jacL[4, 1] += tmp1*(E1dq4_dotL1 + gami*E2dq4_dotL1) + t1*tmp1_dotL1;
  sat_jacL[4, 2] += tmp1*(E1dq4_dotL2 + gami*E2dq4_dotL2) + t1*tmp1_dotL2;
  sat_jacL[4, 3] += tmp1*(E1dq4_dotL3 + gami*E2dq4_dotL3) + t1*tmp1_dotL3;
  sat_jacL[4, 4] += tmp1*(E1dq4_dotL4 + gami*E2dq4_dotL4) + t1*tmp1_dotL4;

  sat_jacR[4, 1] += tmp1*(E1dq4_dotR1 + gami*E2dq4_dotR1) + t1*tmp1_dotR1;
  sat_jacR[4, 2] += tmp1*(E1dq4_dotR2 + gami*E2dq4_dotR2) + t1*tmp1_dotR2;
  sat_jacR[4, 3] += tmp1*(E1dq4_dotR3 + gami*E2dq4_dotR3) + t1*tmp1_dotR3;
  sat_jacR[4, 4] += tmp1*(E1dq4_dotR4 + gami*E2dq4_dotR4) + t1*tmp1_dotR4;


//  println("after second summation")
//  @printit sat_jacR[2, 1]

} // function calcSAT


template <typename Tsol, typename Tmsh, typename Tres>
void calcEulerFlux_diff(Tsol* q, Tres* aux_vars, Tmsh* dir, Tres* Fjac)
{
// calculates the Euler flux in a particular direction at a point
// eqn is the equation type
// q is the vector (of length 4), of the conservative variables at the point
// aux_vars is the vector of auxiliary variables at the point
// dir is a vector of length 2 that specifies the direction
// F is populated with the flux Jacobian
// 2D  only


//  p_dot = params.p_dot
  Tsol p_dot[4];
  press = calcPressure_diff(params, q, p_dot);

  auto fac = 1/q[1];
  auto U = (q[2]*dir[1] + q[3]*dir[2])*fac;
  auto U_dot1 = -(q[2]*dir[1] + q[3]*dir[2])*fac*fac;
  auto U_dot2 = dir[1]*fac;
  auto U_dot3 = dir[2]*fac;

  // F[1] = q[1]*U;
  // F[2] = q[2]*U + dir[1]*press;
  // F[3] = q[3]*U + dir[2]*press;
  // F[4] = (q[4] + press)*U;
  Fjac[1, 1] = U + q[1]*U_dot1;
  Fjac[2, 1] =     q[2]*U_dot1 + dir[1]*p_dot[1];
  Fjac[3, 1] =     q[3]*U_dot1 + dir[2]*p_dot[1];
  Fjac[4, 1] =     q[4]*U_dot1 + press*U_dot1 + U*p_dot[1];

  Fjac[1, 2] =     q[1]*U_dot2;
  Fjac[2, 2] = U + q[2]*U_dot2 + dir[1]*p_dot[2];
  Fjac[3, 2] =     q[3]*U_dot2 + dir[2]*p_dot[2];
  Fjac[4, 2] =     q[4]*U_dot2 + press*U_dot2 + U*p_dot[2];

  Fjac[1, 3] =     q[1]*U_dot3;
  Fjac[2, 3] =     q[2]*U_dot3 + dir[1]*p_dot[3];
  Fjac[3, 3] = U + q[3]*U_dot3 + dir[2]*p_dot[3];
  Fjac[4, 3] =     q[4]*U_dot3 + press*U_dot3 + U*p_dot[3];

  Fjac[1, 4] = 0;
  Fjac[2, 4] = dir[1]*p_dot[4];
  Fjac[3, 4] = dir[2]*p_dot[4];
  Fjac[4, 4] = U + U*p_dot[4];
}


template <typename Tsol>
Tsol calcPressure_diff(Tsol* q, Tsol* p_dot)
{
  // calculate pressure for a node
  // q is a vector of length 4 of the conservative variables

  auto t1 = 1/(q[1]*q[1]);
  auto t2 = q[2]*q[2];
  auto t3 = q[3]*q[3];

  const double gamma_1 = 0.4;
  p_dot[1] = (gamma_1)*( 0.5*(t2*t1 + t3*t1));
  p_dot[2] = -(gamma_1)*(q[2]/q[1]);
  p_dot[3] = -(gamma_1)*(q[3]/q[1]);
  p_dot[4] =   gamma_1;
  return  (gamma_1)*(q[4] - 0.5*(t2 + t3)/q[1]);
}




}

//endif
