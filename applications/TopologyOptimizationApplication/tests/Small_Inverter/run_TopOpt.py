#import kratos core and applications
from numpy import zeros
import numpy as np
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication as ksm
import KratosMultiphysics.TopologyOptimizationApplication as kto
import KratosMultiphysics.LinearSolversApplication as kls
import os
from importlib import import_module
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.TopologyOptimizationApplication import topology_optimizer_factory
from KratosMultiphysics import process_factory
from KratosMultiphysics.KratosUnittest import TestCase
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.StructuralMechanicsApplication import structural_response_function_factory


parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = km.Parameters(parameter_file.read())
parameter_file_pseudo = open("ProjectParametersCompliant.json",'r')
ProjectParameters_pseudo = km.Parameters(parameter_file_pseudo.read())
optimization_file = open("Optimization_Parameters.json",'r')
OptimizationParameters = km.Parameters(optimization_file.read())

#Optimization results to compare for the unittest
results = km.Parameters(""" 
                    {

                    "compliance"                            : 0.0,
                    "number_of_iterations"                  : 0,
                    "volume_fraction"                       : 0.0
                    
                    } 
                    """)

with open("adjoint_max_stress_response_parameters.json",'r') as parameter_sensitivity_stress:
    parameters_sensitivity = km.Parameters(parameter_sensitivity_stress.read())

echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()
current_model = km.Model()
#pseudo
current_model_pseudo = km.Model()
dimension = 2
model_part = current_model.CreateModelPart("Structure")
#pseudo
model_part_pseudo = current_model_pseudo.CreateModelPart("Structure")
solver_module = ProjectParameters["solver_settings"]["solver_type"].GetString()
mod = 'KratosMultiphysics.TopologyOptimizationApplication.topology_optimization_simp_static_solver'

response_function = structural_response_function_factory.CreateResponseFunction("topologyTry", parameters_sensitivity["response_settings"], current_model)
#solver = import_module(mod).CreateSolver(current_model, ProjectParameters["solver_settings"])
#pseudo solver
solver_pseudo= import_module(mod).CreateSolver(current_model_pseudo,ProjectParameters_pseudo["solver_settings"])

model = km.Model()
#with open("adjoint_max_stress_response_parameters.json",'r') as parameter_sensitivity_stress:
 #   parameters_sensitivity = km.Parameters(parameter_sensitivity_stress.read())
#response_function = structural_response_function_factory.CreateResponseFunction("topologyTry", parameters_sensitivity["response_settings"], model)
#model_stress = model.CreateModelPart("Structure")


#solver.AddVariables()
#solver.ImportModelPart()
#solver.AddDofs()
#pseudo
solver_pseudo.AddVariables()
solver_pseudo.ImportModelPart()
solver_pseudo.AddDofs()

#===================================================================

model_part.GetProperties()[1].SetValue(km.YOUNG_MODULUS, 4)
model_part.GetProperties()[1].SetValue(km.POISSON_RATIO, 0.3)
model_part.GetProperties()[1].SetValue(km.DENSITY, 1)
model_part.GetProperties()[1].SetValue(kto.YOUNGS_MODULUS_MIN, 1e-9)
model_part.GetProperties()[1].SetValue(kto.YOUNGS_MODULUS_0, 4)
#pseudo
model_part_pseudo.GetProperties()[1].SetValue(km.YOUNG_MODULUS, 4)
model_part_pseudo.GetProperties()[1].SetValue(km.POISSON_RATIO, 0.3)
model_part_pseudo.GetProperties()[1].SetValue(km.DENSITY, 1)
model_part_pseudo.GetProperties()[1].SetValue(kto.YOUNGS_MODULUS_MIN, 1e-9)
model_part_pseudo.GetProperties()[1].SetValue(kto.YOUNGS_MODULUS_0, 4)

#model_stress.GetProperties()[1].SetValue(km.YOUNG_MODULUS, 3000)
#model_stress.GetProperties()[1].SetValue(km.POISSON_RATIO, 0.3)
#model_stress.GetProperties()[1].SetValue(km.DENSITY, 1)
#model_stress.GetProperties()[1].SetValue(kto.YOUNGS_MODULUS_MIN, 1e-9)
#model_stress.GetProperties()[1].SetValue(kto.YOUNGS_MODULUS_0, 3000)
#===================================================================


cons_law = ksm.LinearElastic3DLaw()
model_part.GetProperties()[1].SetValue(km.CONSTITUTIVE_LAW, cons_law)
model_part_pseudo.GetProperties()[1].SetValue(km.CONSTITUTIVE_LAW, cons_law)
#model_stress.GetProperties()[1].SetValue(km.CONSTITUTIVE_LAW, cons_law)

if(echo_level>1):
    print(model_part)
    for properties in model_part.Properties:
        print(properties)
if(echo_level>1):
    print(model_part_pseudo)
    for properties in model_part_pseudo.Properties:
        print(properties)
        
#if(echo_level>1):
 #   print(model_stress)
  #  for properties in model_stress.Properties:
   #     print(properties)


#obtain the list of the processes to be applied (the process order of execution is important)
#list_of_processes  = process_factory.KratosProcessFactory(current_model).ConstructListOfProcesses( ProjectParameters["processes"]["constraints_process_list"] )
#list_of_processes += process_factory.KratosProcessFactory(current_model).ConstructListOfProcesses( ProjectParameters["processes"]["loads_process_list"] )
#if(ProjectParameters.Has("problem_process_list")):
#    list_of_processes += process_factory.KratosProcessFactory(current_model).ConstructListOfProcesses( ProjectParameters["processes"]["problem_process_list"] )
#if(ProjectParameters.Has("output_process_list")):
 #   list_of_processes += process_factory.KratosProcessFactory(current_model).ConstructListOfProcesses( ProjectParameters["processes"]["output_process_list"] )

#pseudo
list_of_processes_pseudo  = process_factory.KratosProcessFactory(current_model_pseudo).ConstructListOfProcesses( ProjectParameters_pseudo["processes"]["constraints_process_list"] )
list_of_processes_pseudo += process_factory.KratosProcessFactory(current_model_pseudo).ConstructListOfProcesses( ProjectParameters_pseudo["processes"]["loads_process_list"] )
if(ProjectParameters.Has("problem_process_list")):
    list_of_processes_pseudo += process_factory.KratosProcessFactory(current_model_pseudo).ConstructListOfProcesses( ProjectParameters_pseudo["processes"]["problem_process_list"] )
if(ProjectParameters.Has("output_process_list")):
    list_of_processes_pseudo += process_factory.KratosProcessFactory(current_model_pseudo).ConstructListOfProcesses( ProjectParameters_pseudo["processes"]["output_process_list"] )




#===================================================================
           
#if(echo_level>1):
 #   for process in list_of_processes:
  #      print(process)
#for process in list_of_processes:
 #   process.ExecuteInitialize()

#pseudo
for process_pseudo in list_of_processes_pseudo:
    process_pseudo.ExecuteInitialize()
#===================================================================

#computing_model_part = solver.GetComputingModelPart()
#pseudo===========================================================
computing_model_part_pseudo = solver_pseudo.GetComputingModelPart()
#===================================================================
problem_path = os.getcwd()
problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()

# initialize GiD  I/O (gid outputs, file_lists)
#from gid_output_process import GiDOutputProcess
output_settings = ProjectParameters["processes"]["output_configuration"]
#gid_output_process = GiDOutputProcess(computing_model_part,
 #                             problem_name,
  #                            output_settings)

#pseudo
gid_output_process_pseudo = GiDOutputProcess(computing_model_part_pseudo,
                              problem_name,
                              output_settings)

#gid_output_process.ExecuteInitialize()
#solver.Initialize()
#for process in list_of_processes:
 #   process.ExecuteBeforeSolutionLoop()
#gid_output_process.ExecuteBeforeSolutionLoop()

#pseudo=============================================================
gid_output_process_pseudo.ExecuteInitialize()
solver_pseudo.Initialize()
for process_pseudo in list_of_processes_pseudo:
    process_pseudo.ExecuteBeforeSolutionLoop()
gid_output_process_pseudo.ExecuteBeforeSolutionLoop()

response_function.Initialize()
#=====================================================================

def solve_structure(opt_itr, structure_variable):    
    if (structure_variable == 1):
        #for process in list_of_processes:
         #   process.ExecuteInitializeSolutionStep()
        #gid_output_process.ExecuteInitializeSolutionStep()

        #solve problem
        response_function.InitializeSolutionStep()
        response_function.CalculateValue()
        response_function.CalculateGradient()
        stressMax = response_function.GetValue()
        stressNablaE = response_function.GetElementalGradient(ksm.YOUNG_MODULUS_SENSITIVITY)
        ii_3 = 0
        for element_i in model_part.Elements:
            if(ii_3 < OptimizationParameters["optimization_parameters"]["number_of_elements"].GetDouble()):
                Sensitivity_value = stressNablaE[ii_3+1]
                element_i.SetValue(kto.YOUNGS_MODULUS_SENSITIVITY, Sensitivity_value)
                element_i.SetValue(kto.MAX_MEAN_STRESS, stressMax)
                ii_3 += 1
            else:
                ii_3+=1
    
        #for process in list_of_processes:
         #   process.ExecuteFinalizeSolutionStep()
        #gid_output_process.ExecuteFinalizeSolutionStep()
        #for process in list_of_processes:
         #   process.ExecuteFinalizeSolutionStep()
        #for process in list_of_processes:
         #   process.ExecuteBeforeOutputStep()
        #if(gid_output_process.IsOutputStep()):
        #    gid_output_process.PrintOutput()
        #for process in list_of_processes:
         #   process.ExecuteAfterOutputStep() 
    
    if (structure_variable == 0):
        for process_pseudo in list_of_processes_pseudo:
            process_pseudo.ExecuteInitializeSolutionStep()
        gid_output_process_pseudo.ExecuteInitializeSolutionStep() 
        
        #solve problem
        solver_pseudo.InitializeSolutionStep()
        solver_pseudo.SolveSolutionStep()
        solver_pseudo.FinalizeSolutionStep()
    
        for process_pseudo in list_of_processes_pseudo:
            process_pseudo.ExecuteFinalizeSolutionStep()
        gid_output_process_pseudo.ExecuteFinalizeSolutionStep()
        for process_pseudo in list_of_processes_pseudo:
            process_pseudo.ExecuteFinalizeSolutionStep()
        for process_pseudo in list_of_processes_pseudo:
            process_pseudo.ExecuteBeforeOutputStep()
        if(gid_output_process_pseudo.IsOutputStep()):
            gid_output_process_pseudo.PrintOutput()
        for process_pseudo in list_of_processes_pseudo:
            process_pseudo.ExecuteAfterOutputStep()

def FinalizeKSMProcess():
    for process in list_of_processes:
        process.ExecuteFinalize()
    # ending the problem (time integration finished)
    gid_output_process.ExecuteFinalize()
    
def Analyzer(controls, response, opt_itr):
    
    # Save the current number of iterations for comparison reasons of the unittest
    results["number_of_iterations"].SetInt(opt_itr)
    
    # Create object to analyze structure response functions if required
    response_analyzer = kto.StructureResponseFunctionUtilities(model_part)
    linear_solver = km.python_linear_solver_factory.ConstructSolver(ProjectParameters["solver_settings"]["linear_solver_settings"])
    sensitivity_solver = kto.StructureAdjointSensitivityStrategy(model_part, linear_solver,ProjectParameters["solver_settings"]["domain_size"].GetInt())
    
    # Compute objective function value Call the Solid Mechanics Application to compute objective function value
    if(controls["strain_energy"]["calc_func"]):
        # Compute structure solution to get displacement field u
        structure_variable = 1
        solve_structure(opt_itr, structure_variable)
        structure_variable = 0
        solve_structure(opt_itr, structure_variable)

        for element_i in model_part.Elements:
            nn = 0
            disp= [None]*(8*3)
            ID = element_i.Id
            if (ID<OptimizationParameters["optimization_parameters"]["number_of_elements"].GetDouble()):
                for nodes_i in element_i.GetNodes():
                    Id = nodes_i.Id
                    val = model_part_pseudo.Nodes[Id].GetSolutionStepValue(km.DISPLACEMENT)
                    disp[nn*3+0]= val[0]
                    disp[nn*3+1]= val[1]
                    disp[nn*3+2]= val[2]
                    nn += 1
                
                element_i.SetValue(kto.LAMBDA_ADJOINT, disp)

        # Calculate objective function value based on u and save into container
        response["strain_energy"]["func"] =  response_analyzer.ComputeStrainEnergy()

        #Save the compliance of the given iteration as results for the unittest
        results["compliance"].SetDouble(response["strain_energy"]["func"])

    # Compute constraint function value
    if(controls["volume_fraction"]["calc_func"]):
        target_volume_fraction = OptimizationParameters["optimization_parameters"]["initial_volume_fraction"].GetDouble()
        Volume_fraction = response_analyzer.ComputeVolumeFraction()
        response["volume_fraction"]["func"] = Volume_fraction - target_volume_fraction

        #Save the volume fraction of the given iteration as results for the unittest
        results["volume_fraction"].SetDouble(Volume_fraction)

    # Compute sensitivities of objective function
    if(controls["strain_energy"]["calc_grad"]):
        sensitivity_solver.ComputeDisplacementSensitivities()
    # Compute sensitivities of constraint function
    if(controls["volume_fraction"]["calc_grad"]):
        sensitivity_solver.ComputeVolumeFractionSensitivities()
        
    
#    response_function.InitializeSolutionStep()
 #   response_function.CalculateValue()
  #  stressMax = response_function.GetValue()
   # response_function.CalculateGradient()
    #stressMaxNablaShape = response_function.GetNodalGradient(km.SHAPE_SENSITIVITY)
    #stressNablaE = response_function.GetElementalGradient(ksm.YOUNG_MODULUS_SENSITIVITY)
    #response_function.FinalizeSolutionStep()
    
#    print()
 #   print("max stress")
  #  print(np.array(stressMax))
   # print()
    #for i in range(72):
     #   print("Max stress sensitivities w.r.t x, y and z positions for node "+(str(i+1))+":")
      #  print(np.array(stressMaxNablaShape[i+1]))
    #print()

    #for i in range(30):
     #   print("Max stress sensitivities w.r.t E of element "+(str(i+1))+":")
      #  print(np.array(stressNablaE[i+1]))
    #print()

# optimization
optimizer = kto.topology_optimizer_factory.ConstructOptimizer(model_part, model_part_pseudo, OptimizationParameters["optimization_parameters"], Analyzer)
optimizer.optimize()






# Testing the results of the optimization 
TestCase().assertEqual(results["number_of_iterations"].GetInt(), 13)
TestCase().assertAlmostEqual(results["compliance"].GetDouble(), 278.493, 3)
TestCase().assertAlmostEqual(results["volume_fraction"].GetDouble(), 0.50010, 5)

FinalizeKSMProcess()
