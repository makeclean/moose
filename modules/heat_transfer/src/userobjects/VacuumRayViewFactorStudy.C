//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VacuumRayViewFactorStudy.h"
#include "MooseMesh.h"
#include "SideSetPatch.h"

#ifdef MOOSE_XDG_ENABLED
#include "MooseInit.h"
#include "xdg/config.h"
#include "xdg/constants.h"
#include "xdg/error.h"
#include "xdg/mesh_managers.h"
#include "xdg/ray_tracers.h"

#include "libmesh/elem.h"
#include "libmesh/utility.h"

#include <algorithm>
#include <set>
#include <sstream>

namespace
{
/// The sorted global node ids of a side, used as the key of the XDG face map
std::vector<dof_id_type>
sideNodeIds(const libMesh::Elem & elem, const unsigned int side)
{
  // Elem::nodes_on_side returns local node numbers, not global node ids
  std::vector<dof_id_type> node_ids;
  for (const auto local : elem.nodes_on_side(side))
    node_ids.push_back(elem.node_id(local));
  std::sort(node_ids.begin(), node_ids.end());
  return node_ids;
}
}

registerMooseObject("HeatTransferApp", VacuumRayViewFactorStudy);

InputParameters
VacuumRayViewFactorStudy::validParams()
{
  InputParameters params = SideUserObject::validParams();

  MooseEnum qtypes("GAUSS GRID", "GRID");
  params.addParam<MooseEnum>(
      "face_type", qtypes, "The face quadrature rule type used for ray tracing.");

  MooseEnum qorders("CONSTANT FIRST SECOND THIRD FOURTH FIFTH SIXTH SEVENTH EIGHTH NINTH TENTH "
                    "ELEVENTH TWELFTH THIRTEENTH FOURTEENTH FIFTEENTH SIXTEENTH SEVENTEENTH "
                    "EIGHTTEENTH NINTEENTH TWENTIETH",
                    "CONSTANT");
  params.addParam<MooseEnum>(
      "face_order", qorders, "The face quadrature rule order used for ray tracing.");

  params.addParam<unsigned int>(
      "polar_quad_order",
      16,
      "Order of the polar quadrature [polar angle is between ray and normal]. Must be even.");
  params.addParam<unsigned int>(
      "azimuthal_quad_order",
      8,
      "Order of the azimuthal quadrature per quadrant [azimuthal angle is measured in a plane "
      "perpendicular to the normal].");

  params.addParam<std::vector<unsigned int>>(
      "n_patches",
      {},
      "Number of patches per boundary. If empty, every boundary is a single surface. If zero, one "
      "patch is created for every face of the boundary.");

  MultiMooseEnum partitioning(
      "default=-3 metis=-2 parmetis=-1 linear=0 centroid hilbert_sfc morton_sfc");
  partitioning.addValidName("grid");
  params.addParam<MultiMooseEnum>(
      "partitioners",
      partitioning,
      "The partitioner to use for splitting each boundary into patches. Empty (the default) "
      "uses the default partitioner. Only used when 'n_patches' is provided.");

  MultiMooseEnum direction("x y z radial");
  params.addParam<MultiMooseEnum>("centroid_partitioner_directions",
                                  direction,
                                  "The sort direction if using the centroid partitioner.");

  params.addClassDescription(
      "Computes view factors and escape fractions between surfaces radiating into an open "
      "(vacuum) environment by ray tracing with the XDG library.");

  return params;
}

VacuumRayViewFactorStudy::VacuumRayViewFactorStudy(const InputParameters & parameters)
  : SideUserObject(parameters),
    _boundary_ids(_mesh.getBoundaryIDs(getParam<std::vector<BoundaryName>>("boundary"))),
    _n_boundaries(_boundary_ids.size()),
    _n_surfaces(0),
    _fe_face(FEBase::build(_mesh.dimension(), FEType(CONSTANT, MONOMIAL).set_p_refinement(false))),
    _q_face(QBase::build(Moose::stringToEnum<QuadratureType>(getParam<MooseEnum>("face_type")),
                         _mesh.dimension() - 1,
                         Moose::stringToEnum<Order>(getParam<MooseEnum>("face_order")))),
    _aq(std::make_unique<RayTracingAngularQuadrature>(
        3,
        getParam<unsigned int>("polar_quad_order"),
        4 * getParam<unsigned int>("azimuthal_quad_order"),
        0,
        1)),
    _num_dir(_aq->numDirections()),
    _bump(0)
{
  for (const auto i : index_range(_boundary_ids))
    _boundary_index[_boundary_ids[i]] = i;

  // split each participating boundary into patches; a surface is a patch, and the surfaces are
  // ordered boundary by boundary, patch by patch
  std::vector<unsigned int> n_patches = getParam<std::vector<unsigned int>>("n_patches");
  const MultiMooseEnum & partitioners = getParam<MultiMooseEnum>("partitioners");
  const MultiMooseEnum & directions = getParam<MultiMooseEnum>("centroid_partitioner_directions");

  if (n_patches.empty())
    n_patches.assign(_n_boundaries, 1);
  else if (n_patches.size() != _n_boundaries)
    paramError("n_patches", "The number of entries must match the number of boundary entries.");

  if (partitioners.size() != 0 && partitioners.size() != _n_boundaries)
    paramError("partitioners", "The number of entries must match the number of boundary entries.");

  if (directions.size() != 0 && directions.size() != _n_boundaries)
    paramError("centroid_partitioner_directions",
               "The number of entries must match the number of boundary entries.");

  unsigned int offset = 0;
  for (const auto j : index_range(_boundary_ids))
  {
    _patches.push_back(std::make_unique<SideSetPatch>(
        _mesh.getMesh(),
        _boundary_ids[j],
        n_patches[j],
        partitioners.size() ? partitioners[j] : "metis",
        directions.size() ? directions[j] : ""));

    _surface_offsets.push_back(offset);
    _surface_boundary.insert(_surface_boundary.end(), _patches[j]->nPatches(), j);
    offset += _patches[j]->nPatches();
  }
  _n_surfaces = offset;
  _escape_weights.assign(_n_surfaces, 0);
  _vf_weights = ViewFactorMatrix(_n_surfaces);

  if (_mesh.dimension() != 3)
    mooseError("VacuumRayViewFactorStudy only supports three-dimensional meshes.");

  if (!dynamic_cast<const ReplicatedMesh *>(&_mesh.getMesh()))
    mooseError("VacuumRayViewFactorStudy requires a serial (replicated) mesh. "
               "Distributed meshes are not supported.");

  _fe_face->attach_quadrature_rule(_q_face.get());
  _fe_face->get_xyz();
  _fe_face->get_JxW();
  _fe_face->get_normals();
}

void
VacuumRayViewFactorStudy::initialSetup()
{
  // the model is built once on thread 0; the other threaded copies query the
  // primary through viewFactorRow/escapeWeight in finalize
  if (_tid != 0)
    return;

  // hand the whole mesh to XDG: each sideset becomes a surface and the complement
  // of the mesh becomes the void volume the rays travel through
#ifdef XDG_ENABLE_LIBMESH
  // XDG constructs its own libMesh::LibMeshInit unless the host application
  // registers one; a second LibMeshInit would install a second thread-buffered
  // output stream, freeing the one MOOSE installed and leaving libMesh's
  // streams dangling. Register MOOSE's MooseInit before creating any XDG
  // object.
  static bool xdg_registered = false;
  if (!xdg_registered)
  {
    xdg::XDGConfig::config().set_libmesh_external_init(MooseInit::instance());
    xdg_registered = true;
  }
  auto mesh_manager = std::make_shared<xdg::LibMeshManager>(&_mesh.getMesh());
  _xdg = std::make_shared<xdg::XDG>();
  _xdg->set_mesh_manager_interface(mesh_manager);
#ifdef XDG_ENABLE_EMBREE
  _xdg->set_ray_tracing_interface(std::make_shared<xdg::EmbreeRayTracer>());
#else
  mooseError("XDG was built without Embree ray tracing support");
#endif
#else
  mooseError("XDG was built without libMesh support");
#endif
  mesh_manager->init();
  _xdg->prepare_raytracer();

  _void_volume = mesh_manager->implicit_complement();

  // every XDG surface of the model must be one of the participating boundaries,
  // and vice versa: there is no support for non-participating (blocking) surfaces
  // or for sidesets that XDG splits into more than one surface
  std::set<xdg::MeshID> model_surfaces(mesh_manager->surfaces().begin(),
                                       mesh_manager->surfaces().end());
  for (const auto & pair : _boundary_index)
    if (!model_surfaces.count(pair.first))
      mooseError("Boundary ",
                 _mesh.getBoundaryName(pair.first),
                 " does not form a single surface of the XDG model. The vacuum (open geometry) "
                 "view factor calculation does not support sidesets that share interfaces with "
                 "other sidesets.");
  for (const auto surface : model_surfaces)
    if (!_boundary_index.count(surface))
      mooseError("Surface ",
                 surface,
                 " (",
                 _mesh.getBoundaryName(surface),
                 ") of the XDG model does not participate in the radiative exchange. Every "
                 "sideset of the mesh must be listed in the 'boundary' parameter when using "
                 "view_factor_calculator = vacuum_ray_tracing.");

  // map the node ids of every face of the participating surfaces to its XDG face
  // id; this is used to exclude the emitting face from each ray trace
  for (const auto surface : _boundary_ids)
    for (const auto face : mesh_manager->get_surface_faces(surface))
    {
      const auto vertex_ids = mesh_manager->face_vertices(face);
      std::vector<dof_id_type> node_ids(vertex_ids.begin(), vertex_ids.end());
      std::sort(node_ids.begin(), node_ids.end());
      _face_id_by_node_ids[node_ids] = face;
    }

  // map every XDG face to the patch (surface) that contains it, so that a ray hit can be
  // credited to the right patch without looking the face up during every ray trace
  for (const auto & [elem_id, side, sideset] :
       _mesh.getMesh().get_boundary_info().build_active_side_list())
  {
    if (!_boundary_index.count(sideset))
      continue;

    const libMesh::Elem * elem = _mesh.getMesh().elem_ptr(elem_id);
    const auto node_ids = sideNodeIds(*elem, side);

    const auto face_it = _face_id_by_node_ids.find(node_ids);
    if (face_it == _face_id_by_node_ids.end())
      mooseError("A face of boundary ",
                 _mesh.getBoundaryName(sideset),
                 " was not found in the XDG model.");
    _surface_by_face[face_it->second] = surfaceIndex(*elem, side);
  }

  // nudge the ray origin off the emitting surface by the same dilation distance
  // that XDG uses to inflate its bounding boxes, so that the first hit cannot be
  // the emitting face itself
  _bump = std::max(mesh_manager->volume_bounding_box(_void_volume).dilation(), 1e-3);
}

void
VacuumRayViewFactorStudy::initialize()
{
  _vf_weights.clear();
  std::fill(_escape_weights.begin(), _escape_weights.end(), 0);
}

void
VacuumRayViewFactorStudy::execute()
{
  if (_current_boundary_id == Moose::INVALID_BOUNDARY_ID)
    return;

  // each element is computed by the processor that owns it
  if (_current_elem->processor_id() != processor_id())
    return;

  const unsigned int from = surfaceIndex(*_current_elem, _current_side);

  // the XDG model lives on the primary (thread 0) copy of this object
  const VacuumRayViewFactorStudy * model = this;
  if (_tid != 0)
    model = static_cast<const VacuumRayViewFactorStudy *>(primaryThreadCopy());

  // look up the XDG face id of the emitting face so that the source face can be
  // excluded from the ray trace
  const auto node_ids = sideNodeIds(*_current_elem, _current_side);
  const auto face_it = model->_face_id_by_node_ids.find(node_ids);
  mooseAssert(face_it != model->_face_id_by_node_ids.end(),
              "The emitting face is not part of the XDG model.");
  const xdg::MeshID face_id = face_it->second;

  _fe_face->reinit(_current_elem, _current_side);
  const auto & q_point = _fe_face->get_xyz();
  const auto & JxW = _fe_face->get_JxW();
  const auto & normals = _fe_face->get_normals();

  for (const auto qp : index_range(q_point))
  {
    _aq->rotate(normals[qp]);
    const Point origin = q_point[qp] + normals[qp] * model->_bump;

    for (const auto l : make_range(_num_dir))
    {
      const Point & direction = _aq->getDirection(l);
      // the weight of a ray is the outgoing energy carried into its solid angle
      const Real weight = JxW[qp] * (normals[qp] * direction) * _aq->getTotalWeight(l);
      if (weight <= 0)
        continue;

      // a fresh exclude list per ray: ray_fire appends the hit primitive to it
      std::vector<xdg::MeshID> exclude{face_id};
      const xdg::Position pos(origin(0), origin(1), origin(2));
      const xdg::Direction dir(direction(0), direction(1), direction(2));
      const auto [_, hit_surface] = model->_xdg->ray_fire(
          model->_void_volume, pos, dir, xdg::INFTY, xdg::HitOrientation::ANY, &exclude);

      if (hit_surface != xdg::ID_NONE)
        _vf_weights.add(from, libmesh_map_find(model->_surface_by_face, exclude.back()), weight);
      else
        _escape_weights[from] += weight;
    }
  }
}

void
VacuumRayViewFactorStudy::finalize()
{
  // the ray weights are accumulated on the processors that own the elements
  _vf_weights.gatherSum(comm());
  gatherSum(_escape_weights);
}

void
VacuumRayViewFactorStudy::threadJoin(const UserObject & y)
{
  const auto & study = cast_ref<const VacuumRayViewFactorStudy &>(y);
  _vf_weights.merge(study._vf_weights);
  for (const auto i : make_range(_n_surfaces))
    _escape_weights[i] += study._escape_weights[i];
}

unsigned int
VacuumRayViewFactorStudy::surfaceIndex(const libMesh::Elem & elem, const unsigned int side) const
{
  for (const auto j : index_range(_patches))
  {
    const unsigned int patch = _patches[j]->patchIndex(elem.id(), side);
    if (patch != libMesh::invalid_uint)
      return _surface_offsets[j] + patch;
  }

  mooseError("The side ",
             side,
             " of element ",
             elem.id(),
             " is not part of any of the participating boundaries.");
}

unsigned int
VacuumRayViewFactorStudy::surfaceIndex(const xdg::MeshID face) const
{
  return libmesh_map_find(_surface_by_face, face);
}

BoundaryID
VacuumRayViewFactorStudy::surfaceBoundaryID(const unsigned int surface) const
{
  mooseAssert(surface < _n_surfaces, "Surface index out of range.");
  return _boundary_ids[_surface_boundary[surface]];
}

std::string
VacuumRayViewFactorStudy::surfaceName(const unsigned int surface) const
{
  mooseAssert(surface < _n_surfaces, "Surface index out of range.");
  const auto j = _surface_boundary[surface];
  std::stringstream ss;
  ss << _mesh.getBoundaryName(_boundary_ids[j]) << "_" << (surface - _surface_offsets[j]);
  return ss.str();
}

Real
VacuumRayViewFactorStudy::viewFactorWeight(const unsigned int from, const unsigned int to) const
{
  return _vf_weights.get(from, to);
}

const std::vector<ViewFactorMatrix::Entry> &
VacuumRayViewFactorStudy::viewFactorRow(const unsigned int from) const
{
  return _vf_weights.row(from);
}

Real
VacuumRayViewFactorStudy::escapeWeight(const unsigned int from) const
{
  return _escape_weights[from];
}

#endif // MOOSE_XDG_ENABLED
