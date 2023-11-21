import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsTransientThermalTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """
    etalon_value1 = 28.04411163544510063559
    etalon_value2 = 17.55892791313559322
    etalon_value3 = 41.3797035928672316

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass
    
    def test_thermal_heat_flux_2D3N(self):
        test_name = 'test_thermal_heat_flux_2D3N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[37]
        self.assertAlmostEqual(self.etalon_value1, temp)

    def test_thermal_heat_flux_2D6N(self):
        test_name = 'test_thermal_heat_flux_2D6N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[57]
        self.assertAlmostEqual(self.etalon_value1, temp)

    def test_thermal_heat_flux_2D10N(self):
        test_name = 'test_thermal_heat_flux_2D10N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[77]
        self.assertAlmostEqual(self.etalon_value1, temp)

    def test_thermal_heat_flux_2D15N(self):
        test_name = 'test_thermal_heat_flux_2D15N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[97]
        self.assertAlmostEqual(self.etalon_value1, temp)

    def test_thermal_heat_flux_2D4N(self):
        test_name = 'test_thermal_heat_flux_2D4N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[18]
        self.assertAlmostEqual(self.etalon_value2, temp)

    def test_thermal_heat_flux_2D8N(self):
        test_name = 'test_thermal_heat_flux_2D8N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[50]
        self.assertAlmostEqual(self.etalon_value2, temp)

    def test_thermal_heat_flux_2D9N(self):
        test_name = 'test_thermal_heat_flux_2D9N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[63]
        self.assertAlmostEqual(self.etalon_value2, temp)

    def test_thermal_heat_flux_3D4N(self):
        test_name = 'test_thermal_heat_flux_3D4N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[22]
        self.assertAlmostEqual(self.etalon_value3, temp)

    def test_thermal_heat_flux_3D10N(self):
        test_name = 'test_thermal_heat_flux_3D10N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[124]
        self.assertAlmostEqual(self.etalon_value3, temp)
        
    def test_unsteady_thermal_heat_flux_2D6N(self):
        test_name = 'test_unsteady_thermal_heat_flux_2D6N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[57]
        self.assertAlmostEqual(0.46919946397780093, temp)

    def test_unsteady_thermal_heat_flux_2D10N(self):
        test_name = 'test_unsteady_thermal_heat_flux_2D10N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[77]
        self.assertAlmostEqual(0.4674055416030332, temp)

    def test_unsteady_thermal_heat_flux_2D15N(self):
        test_name = 'test_unsteady_thermal_heat_flux_2D15N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[97]
        self.assertAlmostEqual(0.46761403285540487, temp)

    def test_unsteady_thermal_heat_flux_2D4N(self):
        test_name = 'test_unsteady_thermal_heat_flux_2D4N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[18]
        self.assertAlmostEqual(0.12253593527932072, temp)

    def test_unsteady_thermal_heat_flux_2D8N(self):
        test_name = 'test_unsteady_thermal_heat_flux_2D8N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[50]
        self.assertAlmostEqual(0.20716154048406607, temp)

    def test_unsteady_thermal_heat_flux_2D9N(self):
        test_name = 'test_unsteady_thermal_heat_flux_2D9N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[63]
        self.assertAlmostEqual(0.20715104139698065, temp)

    def test_unsteady_thermal_heat_flux_3D4N(self):
        test_name = 'test_unsteady_thermal_heat_flux_3D4N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[22]
        self.assertAlmostEqual(0.8936587648750058, temp)

    def test_unsteady_thermal_heat_flux_3D10N(self):
        test_name = 'test_unsteady_thermal_heat_flux_3D10N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[124]
        self.assertAlmostEqual(1.2294110493528096, temp)
        

if __name__ == '__main__':
    KratosUnittest.main()
