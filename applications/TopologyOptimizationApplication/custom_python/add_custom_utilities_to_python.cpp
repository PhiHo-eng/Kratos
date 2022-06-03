//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//					 Philipp Hofer
//					 Erich Wehrle
//
// ==============================================================================

// External includes 

// Project includes

// Application includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "spaces/ublas_space.h"

// Utilities
#include "custom_utilities/structure_response_function_utilities.h"
#include "custom_utilities/topology_filtering_utilities.h"
#include "custom_utilities/topology_updating_utilities.h"
#include "custom_utilities/topology_extractor_utilities.h"
#include "custom_utilities/topology_smoothing_utilities.h"
#include "custom_utilities/mma_solver.h"
#include "custom_utilities/io_utilities.h"


namespace Kratos
{

namespace Python
{

void AddCustomUtilitiesToPython(pybind11::module& m)

{
    namespace py = pybind11;



    // =============================================================================================================================================
    // Utility Classes
    // =============================================================================================================================================

    py::class_<StructureResponseFunctionUtilities>(m, "StructureResponseFunctionUtilities")
    .def(py::init<ModelPart& >())
    .def("ComputeStrainEnergy", &StructureResponseFunctionUtilities::ComputeStrainEnergy)
    .def("ComputeVolumeFraction", &StructureResponseFunctionUtilities::ComputeVolumeFraction)
    .def("ComputeDisplacementControlledObjective", &StructureResponseFunctionUtilities::ComputeDisplacementControlledObjective)
    .def("ComputeDisplacement", &StructureResponseFunctionUtilities::ComputeDisplacement)
    ;

    py::class_<TopologyFilteringUtilities >(m, "TopologyFilteringUtilities")
    .def(py::init<ModelPart&, const double, const int>())
    .def("ApplyFilterSensitivity", &TopologyFilteringUtilities::ApplyFilterSensitivity)
    .def("ApplyFilterDensity", &TopologyFilteringUtilities::ApplyFilterDensity)
    .def("ApplyFilterStress", &TopologyFilteringUtilities::ApplyFilterStress)
    ;

    py::class_<TopologyUpdatingUtilities >(m, "TopologyUpdatingUtilities")
    .def(py::init<ModelPart&>())
    .def("UpdateDensitiesUsingOCMethod", &TopologyUpdatingUtilities::UpdateDensitiesUsingOCMethod)
    ;
    
    py::class_< MMAAlgorithm >(m, "MMAAlgorithm")
    .def(py::init<ModelPart&, const int, const int>())
    .def("UpdateDensitiesUsingMMAAlgorithm", &MMAAlgorithm::UpdateDensitiesUsingMMAAlgorithm)
    ;

    py::class_<IOUtilities >(m, "IOUtilities" )
    .def(py::init<>())
    .def("SaveOptimizationResults", &IOUtilities::SaveOptimizationResults)
    .def("WriteSurfaceAsSTLFile", &IOUtilities::WriteSurfaceAsSTLFile)
    ;

    py::class_<TopologyExtractorUtilities >(m, "TopologyExtractorUtilities")
    .def(py::init<>())
	.def("ExtractVolumeMesh", &TopologyExtractorUtilities::ExtractVolumeMesh)
	.def("ExtractSurfaceMesh", &TopologyExtractorUtilities::ExtractSurfaceMesh)
	;

	py::class_<TopologySmoothingUtilities >(m, "TopologySmoothingUtilities" )
    .def(py::init<>())
	.def("SmoothMesh", &TopologySmoothingUtilities::SmoothMesh)
	;
}

}  // namespace Python.

} // Namespace Kratos