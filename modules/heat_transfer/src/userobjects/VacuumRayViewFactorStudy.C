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

#ifdef MOOSE_XDG_ENABLED
#include "xdg/constants.h"
#include "xdg/error.h"

#include "libmesh/elem.h"

#include <algorithm>
#include <set>

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

  params.addClassDescription(
      "Computes view factors and escape fractions between surfaces radiating into an open "
      "(vacuum) environment by ray tracing with the XDG library.");

  return params;
}

VacuumRayViewFactorStudy::VacuumRayViewFactorStudy(const InputParameters & parameters)
  : SideUserObject(parameters),
    _boundary_ids(_mesh.getBoundaryIDs(getParam<std::vector<BoundaryName>>("boundary"))),
    _n_boundaries(_boundary_ids.size()),
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
    _vf_weights(_n_boundaries, std::vector<Real>(_n_boundaries, 0)),
    _escape_weights(_n_boundaries, 0)
{
  for (const auto i : index_range(_boundary_ids))
    _boundary_index[_boundary_ids[i]] = i;

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
  // primary through viewFactorWeight/escapeWeight in finalize
  if (_tid != 0)
    return;

  // hand the whole mesh to XDG: each sideset becomes a surface and the complement
  // of the mesh becomes the void volume the rays travel through
  auto mesh_manager = std::make_shared<xdg::LibMeshManager>(&_mesh.getMesh());
  _xdg = std::make_shared<xdg::XDG>();
  _xdg->set_mesh_manager_interface(mesh_manager);
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
                 "view factor calculation does not support sidesets that are split into patches or "
                 "that share interfaces with other sidesets.");
  for (const auto surface : model_surfaces)
    if (!_boundary_index.count(surface))
      mooseError("Surface ",
                 surface,
                 " (",
                 _mesh.getBoundaryName(surface),
                 ") of the XDG model does not participate in the radiative exchange. Every "
                 "sideset of the mesh must be listed in the 'boundary' parameter when using "
                 "view_factor_calculator = vacuum_ray_tracing.");

  for (const auto surface : model_surfaces)
    _surface_to_boundary[surface] = surface;

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

  // nudge the ray origin off the emitting surface by the same dilation distance
  // that XDG uses to inflate its bounding boxes, so that the first hit cannot be
  // the emitting face itself
  _bump = std::max(mesh_manager->volume_bounding_box(_void_volume).dilation(), 1e-3);
}

void
VacuumRayViewFactorStudy::initialize()
{
  for (auto & row : _vf_weights)
    std::fill(row.begin(), row.end(), 0);
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

  const unsigned int from = boundaryIndex(_current_boundary_id);

  // the XDG model lives on the primary (thread 0) copy of this object
  const VacuumRayViewFactorStudy * model = this;
  if (_tid != 0)
    model = static_cast<const VacuumRayViewFactorStudy *>(primaryThreadCopy());

  // look up the XDG face id of the emitting face so that the source face can be
  // excluded from the ray trace
  const auto node_ids_unsigned = _current_elem->nodes_on_side(_current_side);
  std::vector<dof_id_type> node_ids(node_ids_unsigned.begin(), node_ids_unsigned.end());
  std::sort(node_ids.begin(), node_ids.end());
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
        _vf_weights[from]
                   [boundaryIndex(libmesh_map_find(model->_surface_to_boundary, hit_surface))] +=
            weight;
      else
        _escape_weights[from] += weight;
    }
  }
}

void
VacuumRayViewFactorStudy::finalize()
{
  // the ray weights are accumulated on the processors that own the elements
  for (const auto i : make_range(_n_boundaries))
  {
    gatherSum(_vf_weights[i]);
    gatherSum(_escape_weights[i]);
  }
}

void
VacuumRayViewFactorStudy::threadJoin(const UserObject & y)
{
  const auto & study = cast_ref<const VacuumRayViewFactorStudy &>(y);
  for (const auto i : make_range(_n_boundaries))
  {
    _escape_weights[i] += study._escape_weights[i];
    for (const auto j : make_range(_n_boundaries))
      _vf_weights[i][j] += study._vf_weights[i][j];
  }
}

Real
VacuumRayViewFactorStudy::viewFactorWeight(const BoundaryID from_id, const BoundaryID to_id) const
{
  return _vf_weights[boundaryIndex(from_id)][boundaryIndex(to_id)];
}

Real
VacuumRayViewFactorStudy::escapeWeight(const BoundaryID from_id) const
{
  return _escape_weights[boundaryIndex(from_id)];
}

unsigned int
VacuumRayViewFactorStudy::boundaryIndex(BoundaryID id) const
{
  const auto it = _boundary_index.find(id);
  if (it == _boundary_index.end())
    mooseError("Boundary id ",
               id,
               " (",
               _mesh.getBoundaryName(id),
               ") is not listed in the 'boundary' parameter.");
  return it->second;
}

#endif // MOOSE_XDG_ENABLED