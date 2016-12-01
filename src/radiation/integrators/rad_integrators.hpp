#ifndef RADINTEGRATORS_HPP
#define RADINTEGRATORS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file rad_integrators.hpp
//  \brief definitions for RadIntegrator class
//======================================================================================

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../radiation.hpp" // radiation
#ifdef INCLUDE_CHEMISTRY
#include "../../chemistry/species.hpp"
#endif

class MeshBlock;
class ParameterInput;
class Radiation;

//! \class RadIntegrator
//  \brief integrate algorithm for radiative transfer


class RadIntegrator {
  friend class Radiation;
public:
  RadIntegrator(Radiation *prad, ParameterInput *pin);
  ~RadIntegrator();
  
  Radiation *pmy_rad;
  MeshBlock *pmy_mb;

#ifdef INCLUDE_CHEMISTRY
  AthenaArray<Real> col_tot;
  AthenaArray<Real> col_xp, col_xm, col_yp, col_ym, col_zp, col_zm;
  int n_cols_ang; 
  ChemNetwork* pmy_chemnet;
  //calcuate column within each meshblock
  void GetColMB(int direction);
  //calcuate total column and update radiation
  void UpdateRadiation(int direction);
#endif
private:
  Real rad_G0_; //unshielded radiation field strengh, uniform.
};

#endif // RADINTEGRATORS_HPP