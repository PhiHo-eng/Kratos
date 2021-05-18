from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import os

try:
    import numpy as np
    from .MainKratosROM import TestStructuralMechanicsStaticROM
    from .MainKratosHROM import TestStructuralMechanicsStaticHROM
    numpy_available = True
except:
    numpy_available = False

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities

@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class ROMStaticStruct(KratosUnittest.TestCase):
#########################################################################################

    @KratosUnittest.skipUnless(numpy_available, "numpy is required for RomApplication")
    def test_Struct_Static_ROM_2D(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("ProjectParametersROM.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            Simulation = TestStructuralMechanicsStaticROM(model,parameters)
            Simulation.Run()
            ObtainedOutput = Simulation.EvaluateQuantityOfInterest()
            ExpectedOutput = np.load('ExpectedOutputROM.npy')
            NodalArea = Simulation.EvaluateQuantityOfInterest2()

            up=0
            down=0
            for i in range (len(ObtainedOutput)):
                if ExpectedOutput[i] != 0:
                    up += (NodalArea[i]*(   (1  - (ObtainedOutput[i] / ExpectedOutput[i] )    )**2)  )
                    down +=  NodalArea[i]
            l2 = (np.sqrt(up/down)) *100
            self.assertLess(l2, 1e-12)#percent
            # Cleaning
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
            kratos_utilities.DeleteDirectoryIfExisting("vtk_output")
            for file_name in os.listdir(os.getcwd()):
                if file_name.endswith(".bin") or file_name.endswith(".lst") :
                    kratos_utilities.DeleteFileIfExisting(file_name)

    @KratosUnittest.skipUnless(numpy_available, "numpy is required for RomApplication")
    def test_Struct_Static_HROM_2D(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("ProjectParametersHROM.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            Simulation = TestStructuralMechanicsStaticHROM(model,parameters)
            Simulation.Run()
            ObtainedOutput = Simulation.EvaluateQuantityOfInterest()
            ExpectedOutput = np.load('ExpectedOutputHROM.npy')
            NodalArea = Simulation.EvaluateQuantityOfInterest2()
            print(NodalArea)

            UP=0
            DOWN=0
            for i in range (len(ObtainedOutput)):
                if ExpectedOutput[i] != 0:
                    UP += (NodalArea[i]*(   (1  - (ObtainedOutput[i] / ExpectedOutput[i] )    )**2)  )
                    DOWN +=  NodalArea[i]
            L2 = (np.sqrt(UP/DOWN)) *100
            self.assertLess(L2, 1e-12)#percent
            # Cleaning
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
            kratos_utilities.DeleteDirectoryIfExisting("vtk_output")
            for file_name in os.listdir(os.getcwd()):
                if file_name.endswith(".bin") or file_name.endswith(".lst") :
                    kratos_utilities.DeleteFileIfExisting(file_name)


##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
