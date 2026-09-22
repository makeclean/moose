//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ViewFactorBase.h"

#ifdef MOOSE_XDG_ENABLED

// Forward Declarations
class VacuumRayViewFactorStudy;

/**
 * Computes the view factors and escape fractions between surfaces radiating into
 * an open (vacuum) environment. The ray weights of each surface are computed by
 * a VacuumRayViewFactorStudy (ray tracing through XDG); this object normalizes
 * them by the surface area times pi and exposes the escape fractions.
 *
 * View factor normalization is deliberately not supported: radiation that
 * escapes the open geometry must not be redistributed among the participating
 * surfaces.
 */
class VacuumRayViewFactor : public ViewFactorBase
{
public:
  static InputParameters validParams();

  VacuumRayViewFactor(const InputParameters & parameters);

  /// The fraction of radiation leaving this boundary that escapes to the environment
  virtual Real getEscapeFraction(BoundaryID id) const override;

protected:
  /// A purely virtual function called in finalize; fills in the view factors and
  /// escape fractions from the ray tracing study
  virtual void finalizeViewFactor() override;

  /// A purely virtual function called in finalize; joins the thread copies
  virtual void threadJoinViewFactor(const UserObject & y) override;

private:
  /// The ray tracing study that computes the view factor and escape weights
  const VacuumRayViewFactorStudy & _study;

  /// The fraction of radiation that escapes each surface to the environment
  std::vector<Real> _escape_fractions;
};

#endif // MOOSE_XDG_ENABLED