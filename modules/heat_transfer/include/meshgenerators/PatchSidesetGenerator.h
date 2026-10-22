//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MeshGenerator.h"

// Forward declarations

/**
 * Subdivides a sideset into smaller patches each of which is going to be a new sideset
 */
class PatchSidesetGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  PatchSidesetGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

  unsigned int nPatches() const { return _n_patches; }

protected:
  /// returns the names of the _n_patches subdivisions derived from _sideset
  std::vector<BoundaryName> sidesetNameHelper(const std::string & base_name) const;

  std::unique_ptr<MeshBase> & _input;

  /// the number of patches that this sideset generator divides _sideset into
  unsigned int _n_patches;

  /// The sideset that will be subdivided
  const BoundaryName & _sideset_name;

  /// The sideset that will be subdivided
  BoundaryID _sideset;

  /// If true, the original sideset is removed after the patch sidesets are created
  const bool _remove_original_sideset;
};
