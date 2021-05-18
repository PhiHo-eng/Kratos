import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.ConstitutiveModelsApplication
import KratosMultiphysics.DelaunayMeshingApplication
import KratosMultiphysics.UmatApplication
import KratosMultiphysics.PfemApplication
import KratosMultiphysics.ContactMechanicsApplication
import KratosMultiphysics.PfemSolidMechanicsApplication
import MainSolid


model = KratosMultiphysics.Model()
MainSolid.Solution(model).Run()

