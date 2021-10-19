from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *

## This proces sets the value of a scalar variable using the BofangConditionTemperatureProcess.
## In this case, the scalar value is automatically fixed.

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeInputTableNodalYoungModulusProcess(Model, settings["Parameters"])

class ImposeInputTableNodalYoungModulusProcess(Process):

    def __init__(self, Model, settings ):

        Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]
        input_file_name = settings["input_file_name"].GetInt()
        settings.RemoveValue("min_value")
        settings.RemoveValue("max_value")

        self.table = PiecewiseLinearTable()
        with open(input_file_name,'r') as file_name:
            for j, line in enumerate(file_name):
                file_1 = line.split(" ")
                if (len(file_1)) > 1:
                    self.table.AddRow(float(file_1[0]), float(file_1[1]))

        self.process = DamInputTableNodalYoungModulusProcess(model_part, self.table, settings)


    def ExecuteInitialize(self):

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()
