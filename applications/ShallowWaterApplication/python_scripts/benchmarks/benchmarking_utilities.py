from sys import argv
import KratosMultiphysics as KM

from KratosMultiphysics.ShallowWaterApplication.shallow_water_analysis import ShallowWaterAnalysis

class BenchmarkingUtilities:
    def __init__(self, options = ['regular_analysis','convergence_analysis']):
        self.options = options

    @staticmethod
    def RunCase(parameters):
        model = KM.Model()
        analysis = ShallowWaterAnalysis(model, parameters)
        analysis.Run()

    @staticmethod
    def GetParametersFromListOfProcesses(processes_list, name):
        for item in processes_list:
            if item['python_module'].GetString() == name:
                return item['Parameters']

    @staticmethod
    def GetParametersFromListOfModelers(modelers_list, name):
        for item in modelers_list:
            if item['modeler_name'].GetString() == name:
                return item['Parameters']

    @staticmethod
    def ReplaceSettings(parameters, setting, value):
        if isinstance(value, str):
            parameters[setting].SetString(value)
        elif isinstance(value, int):
            parameters[setting].SetInt(value)
        elif isinstance(value, float):
            parameters[setting].SetDouble(value)
        else:
            raise Exception('Unsupported type')

    @classmethod
    def ReplaceProcessSettings(cls, processes_list, python_module, setting, value):
        parameters = cls.GetParametersFromListOfProcesses(processes_list, python_module)
        cls.ReplaceSettings(parameters, setting, value)

    @classmethod
    def ReplaceModelerSettings(cls, processes_list, python_module, setting, value):
        parameters = cls.GetParametersFromListOfModelers(processes_list, python_module)
        cls.ReplaceSettings(parameters, setting, value)

    def ParseArguments(self, argv):
        if len(argv) > 2:
            err_msg  = 'Too many input arguments.\n'
            err_msg += self.Usage()
            raise Exception(err_msg)
        elif len(argv) == 2:
            mode = argv[1]
        else:
            mode = "regular_analysis"
            wrn_msg  = 'Setting the analysis to "{}"\n\n'.format(mode)
            wrn_msg += self.Usage()
            wrn_msg += '\n\n'
            KM.Logger.PrintWarning(wrn_msg)
        return mode

    def Usage(self):
        usage  = 'Usage of this script:\n'
        usage += '> python MainKratos.py "mode"\n'
        usage += 'The possible modes are:\n - '
        usage += '\n - '.join(self.options)
        return usage

    def PrintUsage(self):
        KM.Logger.PrintInfo(self.Usage())

    def PrintUnknownModeMessage(self, mode = None):
        msg  = 'Unknown mode. The specified value is "{}"\n'.format(mode)
        msg += self.Usage()
        KM.Logger.PrintInfo(msg)
