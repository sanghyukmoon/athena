//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file turb.cpp
//! \brief Problem generator for gravitational collapse of turbulent cloud

// C headers

// C++ headers
#include <cmath>
#include <ctime>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../fft/athena_fft.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

#if SELF_GRAVITY_ENABLED != 2
#error "This problem generator requires Multigrid gravity solver."
#endif


int JeansCondition(MeshBlock *pmb);
Real njeans;

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief
//========================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (NON_BAROTROPIC_EOS) {
    std::stringstream msg;
    msg << "This problem generator does not support adiabatic EOS." << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  // In the unit system where [L] = L_J, [M] = M_J, [T] = t_ff,
  // the gravitational constant becomes G = 3*PI/32.
  SetFourPiG(4.*PI*(3.*PI/32.));

  // turb_flag is initialzed in the Mesh constructor to 0 by default;
  // turb_flag = 1 for decaying turbulence
  // turb_flag = 2 for impulsively driven turbulence
  // turb_flag = 3 for continuously driven turbulence
  turb_flag = pin->GetInteger("turbulence","turb_flag");
  if (turb_flag != 0) {
#ifndef FFT
    std::stringstream msg;
    msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << std::endl
        << "non zero Turbulence flag is set without FFT!" << std::endl;
    ATHENA_ERROR(msg);
    return;
#endif
  }

  if (adaptive) {
    njeans = pin->GetReal("mesh", "njeans");
    EnrollUserRefinementCondition(JeansCondition);
  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IDN,k,j,i) = 1.0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
      }
    }
  }
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
}

// Jeans Condition
int JeansCondition(MeshBlock *pmb) {
  Real njmin = 1e300;
  const Real dx = pmb->pcoord->dx1f(0); // assuming uniform cubic cells
  const Real cs = pmb->peos->GetIsoSoundSpeed();
  const Real fac = std::sqrt(PI) / dx;
  for (int k = pmb->ks; k<=pmb->ke; ++k) {
    for (int j = pmb->js; j<=pmb->je; ++j) {
      for (int i = pmb->is; i<=pmb->ie; ++i) {
        Real nj = fac * cs / std::sqrt(pmb->phydro->w(IDN,k,j,i));
        njmin = std::min(njmin, nj);
      }
    }
  }
  if (njmin < njeans)
    return 1;
  if (njmin > njeans * 2.5)
    return -1;
  return 0;
}
