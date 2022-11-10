//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file isothermal_hydro.cpp
//! \brief implements functions in class EquationOfState for isothermal hydrodynamics`

// C headers

// C++ headers
#include <cmath>   // sqrt()
#include <string>     // std::string

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "eos.hpp"

// EquationOfState constructor

EquationOfState::EquationOfState(MeshBlock *pmb, ParameterInput *pin) :
    pmy_block_(pmb),
    iso_sound_speed_{pin->GetReal("hydro", "iso_sound_speed")},  // error if missing!
    density_floor_{pin->GetOrAddReal("hydro", "dfloor", std::sqrt(1024*float_min) )},
    scalar_floor_{pin->GetOrAddReal("hydro", "sfloor", std::sqrt(1024*float_min))},
    floor_method_{pin->GetOrAddString("hydro", "floor_method", "floor")} {}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
//!           const AthenaArray<Real> &prim_old, const FaceField &b,
//!           AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Coordinates *pco,
//!           int il, int iu, int jl, int ju, int kl, int ku)
void EquationOfState::ConservedToPrimitive(
    AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old, const FaceField &b,
    AthenaArray<Real> &prim, AthenaArray<Real> &bcc,
    Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku) {
  unsigned int nbad_d = 0;
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real& u_d  = cons(IDN,k,j,i);
        Real& u_m1 = cons(IM1,k,j,i);
        Real& u_m2 = cons(IM2,k,j,i);
        Real& u_m3 = cons(IM3,k,j,i);

        Real& w_d  = prim(IDN,k,j,i);
        Real& w_vx = prim(IVX,k,j,i);
        Real& w_vy = prim(IVY,k,j,i);
        Real& w_vz = prim(IVZ,k,j,i);

        // apply density floor, without changing momentum or energy
        nbad_d += CorrectBadCells(cons,k,j,i,il,iu,jl,ju,kl,ku,floor_method_);
        w_d = u_d;

        Real di = 1.0/u_d;
        w_vx = u_m1*di;
        w_vy = u_m2*di;
        w_vz = u_m3*di;
      }
    }
  }
  if (nbad_d > 0)
    std::cerr << "[Bad cell correction] rank " << Globals::my_rank
      << ": " << nbad_d << " bad cells are corrected." << std::endl;

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
//!           const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
//!           int il, int iu, int jl, int ju, int kl, int ku);
//! \brief Converts primitive variables into conservative variables

void EquationOfState::PrimitiveToConserved(
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bc,
    AthenaArray<Real> &cons, Coordinates *pco,
    int il, int iu, int jl, int ju, int kl, int ku) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real& u_d  = cons(IDN,k,j,i);
        Real& u_m1 = cons(IM1,k,j,i);
        Real& u_m2 = cons(IM2,k,j,i);
        Real& u_m3 = cons(IM3,k,j,i);

        const Real& w_d  = prim(IDN,k,j,i);
        const Real& w_vx = prim(IVX,k,j,i);
        const Real& w_vy = prim(IVY,k,j,i);
        const Real& w_vz = prim(IVZ,k,j,i);

        u_d = w_d;
        u_m1 = w_vx*w_d;
        u_m2 = w_vy*w_d;
        u_m3 = w_vz*w_d;
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::SoundSpeed(Real dummy_arg[NHYDRO])
//! \brief returns isothermal sound speed

Real EquationOfState::SoundSpeed(const Real dummy_arg[NHYDRO]) {
  return iso_sound_speed_;
}

//---------------------------------------------------------------------------------------
//! \fn void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim, int k, int j,
//!                                                 int i)
//! \brief Apply density floor to reconstructed L/R cell interface states

void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim, int k, int j, int i) {
  Real& w_d  = prim(IDN,i);

  // apply density floor
  w_d = (w_d > density_floor_) ?  w_d : density_floor_;

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::ApplyPrimitiveConservedFloors(AthenaArray<Real> &prim,
//!           AthenaArray<Real> &cons, FaceField &b, int k, int j, int i) {
//! \brief Apply pressure (prim) floor and correct energy (cons) (typically after W(U))

void EquationOfState::ApplyPrimitiveConservedFloors(
    AthenaArray<Real> &prim, AthenaArray<Real> &cons, AthenaArray<Real> &bcc,
    int k, int j, int i) {
  Real& w_d  = prim(IDN,k,j,i);
  Real& u_d  = cons(IDN,k,j,i);

  // apply (prim) density floor, without changing momentum or energy
  w_d = (w_d > density_floor_) ?  w_d : density_floor_;
  // ensure cons density matches
  u_d = w_d;

  return;
}

//---------------------------------------------------------------------------------------
//! \fn void EquationOfState::CorrectBadCells(AthenaArray<Real> &cons, int k, int j,
//!           int i, int il, int ku, int jl, int ju, int kl, int ku, std::string method)
//! \brief Apply density floor after integration at each stage

int EquationOfState::CorrectBadCells(AthenaArray<Real> &cons, int k, int j, int i,
    int il, int iu, int jl, int ju, int kl, int ku, std::string method) {
  Real& u_d = cons(IDN,k,j,i);
  Real u_d_old = cons(IDN,k,j,i);

  if (method == "floor") {
    // apply (cons) density floor, without changing momentum or energy
    if (u_d <= density_floor_) {
      u_d = density_floor_;
      std::cerr << "[Bad cell correction] rank " << Globals::my_rank
        << ": density floor applied. old=" << u_d_old << " new=" << u_d << std::endl;
      return 1;
    }
  } else if (method == "neighbor") {
    // apply neighbor averaging
    // (SM) WARNING: The loop calling this function will not be vectorized
    if (u_d <= density_floor_) {
      Real u_d_avg = 0;
      int koff[] = {1,-1,0,0,0,0};
      int joff[] = {0,0,1,-1,0,0};
      int ioff[] = {0,0,0,0,1,-1};
      int n_neighbors = 0;
      for (int idx=0; idx<6; ++idx) {
        int kk = k + koff[idx];
        int jj = j + joff[idx];
        int ii = i + ioff[idx];
        // skip indices outside the range
        if ((ii<il) || (ii>iu) || (jj<jl) || (jj>ju) || (kk<kl) || (kk >ku))
          continue;
        n_neighbors++;
        u_d_avg += cons(IDN,kk,jj,ii);
      }
      u_d = u_d_avg / (Real)n_neighbors;
      std::cerr << "[Bad cell correction] rank " << Globals::my_rank
        << ": density floor applied. old=" << u_d_old << " new=" << u_d << std::endl;
      return 1;
    }
  }
  return 0;
}
