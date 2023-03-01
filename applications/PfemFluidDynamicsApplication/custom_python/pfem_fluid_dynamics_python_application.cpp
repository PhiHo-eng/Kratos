//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes
#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_modelers_to_python.h"
#include "custom_python/add_custom_bounding_to_python.h"

#include "includes/define.h"
#include "includes/define_python.h"

#include "pfem_fluid_dynamics_application.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

PYBIND11_MODULE(KratosPfemFluidDynamicsApplication, m)
{

  py::class_<KratosPfemFluidDynamicsApplication,
             KratosPfemFluidDynamicsApplication::Pointer,
             KratosApplication>(m, "KratosPfemFluidDynamicsApplication")
      .def(py::init<>());

  AddCustomProcessesToPython(m);
  AddCustomUtilitiesToPython(m);
  AddCustomStrategiesToPython(m);
  AddCustomConstitutiveLawsToPython(m);
  AddCustomModelersToPython(m);
  AddCustomBoundingToPython(m);

  //registering variables in python ( if must to be seen from python )

  // some post process variables + stress invariants
  // KRATOS_REGISTER_IN_PYTHON_VARIABLE( M_MODULUS )
  // KRATOS_REGISTER_IN_PYTHON_VARIABLE(PATCH_INDEX);
  // KRATOS_REGISTER_IN_PYTHON_VARIABLE(NORMVELOCITY);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NO_MESH);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EULERIAN_INLET);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGIAN_INLET);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FREESURFACE);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PREVIOUS_FREESURFACE);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INITIAL_DELTA_TIME);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CURRENT_DELTA_TIME);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TIME_INTERVAL_CHANGED);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ISOLATED_NODE);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_H_WALL);

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MECHANICAL_DISSIPATION);

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, YIELDED);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FLOW_INDEX);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, YIELD_SHEAR);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ADAPTIVE_EXPONENT);

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COHESION);

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, STATIC_FRICTION);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DYNAMIC_FRICTION);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INERTIAL_NUMBER_ZERO);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, GRAIN_DIAMETER);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, GRAIN_DENSITY);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, REGULARIZATION_COEFFICIENT);

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESSURE_VELOCITY);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESSURE_ACCELERATION);

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_ERROR_XX);

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_CAUCHY_STRESS);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_DEVIATORIC_CAUCHY_STRESS);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_SFD_NEIGHBOURS);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_SFD_NEIGHBOURS_ORDER);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_DEFORMATION_GRAD);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_DEFORMATION_GRAD_VEL);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_SPATIAL_DEF_RATE);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_VOLUMETRIC_DEF_RATE);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_EQUIVALENT_STRAIN_RATE);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_MEAN_MESH_SIZE);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_TAU);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_FREESURFACE_AREA);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VOLUMETRIC_COEFFICIENT);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEVIATORIC_COEFFICIENT);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INTERFACE_NODE);

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_NODAL_VOLUME);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_NODAL_CAUCHY_STRESS);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_NODAL_SFD_NEIGHBOURS);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_NODAL_SFD_NEIGHBOURS_ORDER);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_NODAL_DEFORMATION_GRAD);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_NODAL_DEFORMATION_GRAD_VEL);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_NODAL_SPATIAL_DEF_RATE);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_NODAL_VOLUMETRIC_DEF_RATE);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_NODAL_EQUIVALENT_STRAIN_RATE);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_NODAL_MEAN_MESH_SIZE);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_DENSITY);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_NODAL_TAU);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_NODAL_FREESURFACE_AREA);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_VOLUMETRIC_COEFFICIENT);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_DEVIATORIC_COEFFICIENT);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_INTERFACE_NODE);
}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
