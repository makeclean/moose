//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SideUserObject.h"
#include "MooseEnum.h"

#include <petscmat.h>
#include <petscksp.h>

// Forward Declarations
class Function;
class ViewFactorMatrix;

/**
 * GrayLambertSurfaceRadiationBase computes the heat flux on a set of surfaces
 * in radiative heat transfer with each other.
 *
 * The surfaces are the patches of the participating boundaries. The patches are either materialized
 * as sidesets (in which case every sideset is one surface) or indexed internally, in which case the
 * boundaries listed in 'boundary' are split into 'n_patches' patches each.
 */
class GrayLambertSurfaceRadiationBase : public SideUserObject
{
public:
  static InputParameters validParams();

  GrayLambertSurfaceRadiationBase(const InputParameters & parameters);
  virtual ~GrayLambertSurfaceRadiationBase();

  virtual void execute() override;
  virtual void initialize() override;
  virtual void finalize() override;
  bool checkVariableBoundaryIntegrity() const override { return false; }

  /// Define enum for boundary type
  enum RAD_BND_TYPE
  {
    VARIABLE_TEMPERATURE = 0,
    FIXED_TEMPERATURE = 4,
    ADIABATIC = 8
  };

  ///@{ public interface of this UserObject; boundary-keyed quantities are area-weighted
  ///    averages over the patches of the boundary when the patches are internal
  Real getSurfaceIrradiation(BoundaryID id) const;
  Real getSurfaceHeatFluxDensity(BoundaryID id) const;
  Real getSurfaceTemperature(BoundaryID id) const;
  Real getSurfaceRadiosity(BoundaryID id) const;
  Real getSurfaceEmissivity(BoundaryID id) const;
  Real getViewFactor(BoundaryID from_id, BoundaryID to_id) const;
  std::set<BoundaryID> getSurfaceIDs() const;
  ///@}

  ///@{ per-patch quantities, for applying a flux on the side that is being integrated
  unsigned int surfaceIndex(const libMesh::Elem & elem, unsigned int side) const;
  Real getSurfaceIrradiation(const libMesh::Elem & elem, unsigned int side) const;
  Real getSurfaceHeatFluxDensity(const libMesh::Elem & elem, unsigned int side) const;
  Real getSurfaceEmissivity(const libMesh::Elem & elem, unsigned int side) const;
  ///@}

protected:
  virtual void threadJoin(const UserObject & y) override;

  /// a purely virtual function that defines where view factors come from
  virtual const ViewFactorMatrix * setViewFactors() = 0;

  /// returns the fraction of radiation that escapes each surface to the environment
  virtual std::vector<Real> setEscapeFractions();

  /// Gets the index within \c _side_id_index of the given boundary ID
  unsigned int getSideIDIndex(BoundaryID id) const;

  /// Gets the parameter index of the given boundary ID
  unsigned int getBoundaryIndex(BoundaryID id) const;

  /// The surface index range [begin, end) of the boundary with parameter index \p j
  std::pair<unsigned int, unsigned int> surfaceRange(unsigned int j) const;

  /// Stefan-Boltzmann constant
  const Real _sigma_stefan_boltzmann;

  /// number of active boundary ids
  unsigned int _n_sides;

  /// the coupled temperature variable
  const VariableValue & _temperature;

  /// the participating boundary ids, in the order of the 'boundary' parameter
  std::vector<BoundaryID> _boundary_ids;

  /// map from boundary id to the index in the 'boundary' parameter
  std::map<BoundaryID, unsigned int> _boundary_index;

  /// whether the participating boundaries are split into patches internally
  bool _internal_patches;

  /// the number of patches of each participating boundary
  std::vector<unsigned int> _n_patches_per_boundary;

  /// the first surface index of each participating boundary
  std::vector<unsigned int> _surface_offsets;

  /// the index in the 'boundary' parameter of each surface
  std::vector<unsigned int> _surface_boundary;

  /// (element id, side) to surface index map, only filled for internal patches
  std::map<std::pair<dof_id_type, unsigned int>, unsigned int> _side_surface;

  /// the fixed temperature function of each surface, null if the surface is not isothermal
  std::vector<const Function *> _fixed_side_temperature;

  /// the radiosity of each surface
  std::vector<Real> _radiosity;

  /// the heat flux density qdot
  std::vector<Real> _heat_flux_density;

  /// the average temperature: this could be important for adiabatic walls
  std::vector<Real> _side_temperature;

  /// the type of the side, allows lookup index -> type
  std::vector<enum RAD_BND_TYPE> _side_type;

  /// side id to index map, side ids can have holes or be out of order
  std::map<BoundaryID, unsigned int> _side_id_index;

  /// the area by participating side set
  std::vector<Real> _areas;

  /// the average value of sigma * eps * T^4
  std::vector<Real> _beta;

  /// the irradiation into each surface
  std::vector<Real> _surface_irradiation;

  /// the fraction of radiation that escapes from each surface to the environment
  std::vector<Real> _escape_fractions;

  /// the ambient environment of the surfaces
  const MooseEnum _environment;

  /// the temperature of the environment if it is a black body
  const Function * _environment_temperature;

  /// constant emissivity for each surface
  std::vector<const Function *> _emissivity;

  /// the set of adiabatic boundaries
  std::set<unsigned int> _adiabatic_side_ids;

  /// the view factors which are set by setViewFactors by derived classes
  const ViewFactorMatrix * _view_factors;

  /// matrix A = I - F o diag(c) of the radiosity system, created in initialize()
  Mat _mat;

  /// GMRES solver for the radiosity system
  KSP _ksp;

  /// ILU(0) preconditioner of the radiosity KSP; owned by the KSP
  PC _pc;

  /// right-hand side of the radiosity system
  Vec _rhs;

  /// solution of the radiosity system
  Vec _radiosity_vec;
};
