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
#include "ViewFactorMatrix.h"

// Forward Declarations

/**
 * A base class for automatic computation of view factors between sidesets
 */
class ViewFactorBase : public SideUserObject
{
public:
  static InputParameters validParams();

  ViewFactorBase(const InputParameters & parameters);

  /// How to normalize view factors
  enum class NormalizationMethod
  {
    LAGRANGE_MULTIPLIER = 0,
    INVERSE_ROW_SUM = 1
  };

  ///@{ public interface for obtaining view factors
  Real getViewFactor(BoundaryID from_id, BoundaryID to_id) const;
  Real getViewFactor(BoundaryName from_name, BoundaryName to_name) const;
  /// The full view factor matrix, sparse, indexed by the boundary parameter order
  const ViewFactorMatrix & viewFactorMatrix() const;
  ///@}

  /// The number of surfaces (patches) the view factors are defined on. This is one per
  /// boundary unless the derived object splits the boundaries into patches internally.
  virtual unsigned int nSurfaces() const { return _n_sides; }

  /// Whether the surfaces are patches that are indexed internally instead of separate
  /// boundary ids
  virtual bool usesInternalPatches() const { return false; }

  /// The index of the surface that contains the given mesh side
  virtual unsigned int surfaceIndex(const libMesh::Elem & elem, unsigned int side) const;

  /// The fraction of radiation leaving this boundary that escapes to the environment.
  /// The default is zero for closed geometries; open geometries override this.
  virtual Real getEscapeFraction(BoundaryID /*id*/) const { return 0; }

  /// The escape fraction of a surface, by surface index
  virtual Real getEscapeFraction(unsigned int /*surface*/) const { return 0; }

  /// A name for a surface, used in reporting. Defaults to the boundary name.
  virtual std::string surfaceName(unsigned int surface) const;

  /// The boundary id a surface belongs to, used in reporting
  virtual BoundaryID surfaceBoundaryID(unsigned int surface) const;

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override final;

protected:
  virtual void threadJoin(const UserObject & y) override final;

  /// this function computes the deviation from reciprocity defined as Fij - Aj/Ai * Fji
  Real devReciprocity(unsigned int i, unsigned int j) const;

  /// this function computes the maximum absolute value of the deviation from reciprocity
  Real maxDevReciprocity() const;

  /// sum of a row in the view factor matrix
  Real viewFactorRowSum(unsigned int i) const;

  /// maximum deviation of any view factor row sum from 1
  Real maxDevRowSum() const;

  /// normalizes view factors by Lagrange multiplier method
  void normalizeUsingLagrangeMultiplier();

  /// normalizes view factors by inverse row sum
  void normalizeUsingInverseRowSum();

  /// helper for finding index of correction for i,j-th entry
  unsigned int indexHelper(unsigned int i, unsigned int j) const;

  /// helper function to get the index from the boundary name
  unsigned int getSideNameIndex(std::string name) const;

  /// this function checks & normalizes view factors to sum to one, this is not always
  void checkAndNormalizeViewFactor();

  /// a purely virtural function called in finalize, must be overriden by derived class
  virtual void finalizeViewFactor() = 0;

  /// a purely virtual function called in finalize, must be overriden by derived class
  virtual void threadJoinViewFactor(const UserObject & y) = 0;

  /// number of boundaries of this side uo
  unsigned int _n_sides;

  /// area of the sides i
  std::vector<Real> _areas;

  /// view factor tolerance
  const Real _view_factor_tol;

  /// whether to normalize view factors so vf[from][:] sums to one
  const bool _normalize_view_factor;

  /// How to normalize view factors
  const NormalizationMethod _normalization_method;

  // whether to print view factor information
  const bool _print_view_factor_info;

  /// the view factor from side i to side j, sparse
  ViewFactorMatrix _view_factors;

  /// boundary name to index map
  std::unordered_map<std::string, unsigned int> _side_name_index;
};
