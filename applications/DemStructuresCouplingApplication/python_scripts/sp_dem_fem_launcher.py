import sys
import KratosMultiphysics as Kratos
from Kratos import Logger
import KratosMultiphysics.DEMApplication as Dem
import KratosMultiphysics.ExternalSolversApplication as ExternalSolvers
import KratosMultiphysics.StructuralMechanicsApplication as Structural
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem

from KratosMultiphysics.DemStructuresCouplingApplication.sp_dem_fem_coupling_algorithm import SPAlgorithm

SPAlgorithm().Run()
