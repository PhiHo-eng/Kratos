import KratosMultiphysics.sympy_fe_utilities as KratosSympy

## Computation of the Source Matrix
def ComputeSourceMatrix(dofs, mass_source, body_force, heat_source, params):
    """This function calculates the source matrix"""

    print("\nCompute Source Matrix \n")
    dim = params.dim	# Spatial dimensions

    ## Source fields definition
    rho = dofs[0]           # Density
    m = mass_source         # Mass source (kg/m³s)
    f = body_force.copy()	# Body force vector
    r = heat_source			# Heat source/sink

    '''
    S - Reactive matrix definition
     m/rho 0  0  0
     fx    0  0  0
     fy    0  0  0
     r     fx fy 0
    '''
    S = KratosSympy.zeros(dim+2, dim+2) # Reactive matrix (source terms)
    S[0, 0] = m / rho # Mass source term (divided by rho as S matrix will be later on multiplied by U vector)
    for i in range(1, dim+1): # Body force
        S[i, 0] = f[i-1]
        S[dim+1, i] = f[i-1]
    S[dim+1, 0] = r # Heat source term

    return S

def PrintSourceMatrix(S,params):
    """Auxiliary function to print the source matrix (S)"""

    dim = params.dim	# Spatial dimensions
    print("The source term matrix is:\n")
    for i in range (0,dim+2):
        for j in range (0,dim+2):
            print("S[",i,",",j,"]=",S[i,j],"\n")

    return 0
