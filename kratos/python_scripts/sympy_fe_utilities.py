import sympy
import re

def DefineMatrix(name, m, n):
    """ This method defines a symbolic matrix

    Keyword arguments:
    name -- Name of variables.
    m -- Number of rows.
    n -- Number of columns.
    """
    return sympy.Matrix(m, n, lambda i, j: sympy.var("{name}_{i}_{j}".format(name=name, i=i, j=j)))

def DefineSymmetricMatrix(name, m, n):
    """ This method defines a symbolic symmetric matrix

    Keyword arguments:
    name -- Name of variables.
    m -- Number of rows.
    n -- Number of columns.
    """
    return sympy.Matrix(m, n, lambda i, j:
        sympy.var("{name}_{i}_{j}".format(name=name, i=min(i,j), j=max(i,j))))

def DefineVector( name, m):
    """ This method defines a symbolic vector

    Keyword arguments:
    name -- Name of variables.
    m -- Number of components.
    """
    return sympy.Matrix(m, 1, lambda i,_: sympy.var("{name}_{i}".format(name=name, i=i)))

def DefineShapeFunctions(nnodes, dim, impose_partion_of_unity=False):
    """ This method defines shape functions and derivatives
    Note that partition of unity is imposed
    the name HAS TO BE --> N and DN

    Keyword arguments:
    nnodes -- Number of nodes
    dim -- Dimension of the space
    impose_partion_of_unity -- Impose the partition of unity

    Note that partition of unity is imposed the name HAS TO BE --> N and DN
    """
    DN = DefineMatrix('DN', nnodes, dim)
    N = DefineVector('N', nnodes)

    #impose partition of unity
    if impose_partion_of_unity:
        N[nnodes-1] = 1
        for i in range(nnodes-1):
            N[nnodes-1] -= N[i]

        DN[nnodes-1,:] = -DN[0,:]
        for i in range(1,nnodes-1):
            DN[nnodes-1,:] -= DN[i,:]

    return N, DN

def StrainToVoigt(M):
    """ This method transform the strains matrix to Voigt notation

    Keyword arguments:
    M -- The strain matrix
    """
    if(M.shape[0] == 2):
        vm = sympy.Matrix(3, 1, lambda _: 0.0)
        vm[0,0] = M[0,0]
        vm[1,0] = M[1,1]
        vm[2,0] = 2.0*M[0,1]
    elif(M.shape[0] == 3):
        raise NotImplementedError()
    return vm

def MatrixB(DN):
    """ This method defines the deformation matrix B

    Keyword arguments:
    DN -- The shape function derivatives
    """
    dim = DN.shape[1]
    if dim == 2:
        strain_size = 3
        nnodes = DN.shape[0]
        B = sympy.Matrix(sympy.zeros(strain_size, nnodes*dim))
        for i in range(nnodes):
            for _ in range(dim):
                B[0, i*dim] = DN[i,0]
                B[0, i*dim+1] = 0
                B[1, i*dim] = 0
                B[1, i*dim+1] = DN[i,1]
                B[2, i*dim] = DN[i,1]
                B[2, i*dim+1] = DN[i,0]
    elif dim == 3:
        strain_size = 6
        nnodes = DN.shape[0]
        B = sympy.Matrix(sympy.zeros(strain_size, nnodes*dim))
        for i in range(nnodes):
            B[0, i*3 ] = DN[i, 0]
            B[1, i*3 + 1] = DN[i, 1]
            B[2, i*3 + 2] = DN[i, 2]
            B[3, i*3] = DN[i, 1]
            B[3, i*3 + 1] = DN[i, 0]
            B[4, i*3 + 1] = DN[i, 2]
            B[4, i*3 + 2] = DN[i, 1]
            B[5, i*3] = DN[i, 2]
            B[5, i*3 + 2] = DN[i, 0]
    else:
        print("dimension asked in Matrix B is ",dim)
        raise ValueError("wrong dimension")
    return B

def grad_sym_voigtform(DN, x):
    """ This method defines a symmetric gradient

    Keyword arguments:
    DN -- The shape function derivatives
    x -- The variable to compute the gradient
    """
    [nnodes, dim] = x.shape

    B = MatrixB(DN)

    # Put the x components one after the other in a vector
    xvec = sympy.Matrix(sympy.zeros(B.shape[1], 1))
    for i in range(nnodes):
        for k in range(dim):
            xvec[i*dim+k] = x[i,k]

    return sympy.simplify(B*xvec)

def DfjDxi(DN,f):
    """ This method defines a gradient. This returns a matrix D such that D(i,j) = D(fj)/D(xi)

    Keyword arguments:
    DN -- The shape function derivatives
    f-- The variable to compute the gradient
    """
    return sympy.simplify(DN.transpose()*f)

def DfiDxj(DN,f):
    """ This method defines a gradient This returns a matrix D such that D(i,j) = D(fi)/D(xj)

    Keyword arguments:
    DN -- The shape function derivatives
    f -- The variable to compute the gradient
    """
    return (DfjDxi(DN,f)).transpose()

def div(DN,x):
    """ This method defines the divergence

    Keyword arguments:
    DN -- The shape function derivatives
    x -- The variable to compute the gradient
    """
    if(DN.shape != x.shape):
        raise Exception("shapes are not compatible")

    div_x = 0
    for i in range(DN.shape[0]):
        for k in range(DN.shape[1]):
            div_x += DN[i,k]*x[i,k]

    return sympy.Matrix([sympy.simplify(div_x)])

def SubstituteMatrixValue(where_to_substitute, what_to_substitute, substituted_value):
    """ This method substitutes values into a matrix

    Keyword arguments:
    where_to_substitute -- Coordinates where to substitute
    what_to_substitute -- Components to substitute
    substituted_value -- Variable to substitute
    """
    for lll in range(where_to_substitute.shape[0]):
        for kkk in range(where_to_substitute.shape[1]):
            tmp = where_to_substitute[lll, kkk]
            for i in range(what_to_substitute.shape[0]):
                for j in range(what_to_substitute.shape[1]):
                    tmp = tmp.subs(what_to_substitute[i,j], substituted_value[i,j])

            where_to_substitute[lll, kkk] = tmp

    return where_to_substitute

def SubstituteScalarValue(where_to_substitute, what_to_substitute, substituted_value):
    """ This method substitutes values into a scalar

    Keyword arguments:
    where_to_substitute -- Coordinates where to substitute
    what_to_substitute -- Components to substitute
    substituted_value -- Variable to substitute
    """
    for lll in range(where_to_substitute.shape[0]):
        tmp  = where_to_substitute[lll]
        tmp = tmp.subs( what_to_substitute, substituted_value)
        where_to_substitute[lll] = tmp
    return where_to_substitute

def Compute_RHS(functional, testfunc, do_simplifications=False):
    """ This computes the RHS vector

    Keyword arguments:
    functional -- The functional to derivate
    testfunc -- The test functions
    do_simplifications -- If apply simplifications
    """
    rhs = sympy.Matrix(sympy.zeros(testfunc.shape[0],1))
    for i in range(testfunc.shape[0]):
        rhs[i] = sympy.diff(functional[0,0], testfunc[i])

        if do_simplifications:
            rhs[i] = sympy.simplify(rhs[i])

    return rhs

def Compute_LHS(rhs, testfunc, dofs, do_simplifications=False):
    """ This computes the LHS matrix

    Keyword arguments:
    rhs -- The RHS vector
    testfunc -- The test functions
    dofs -- The dofs vectors
    do_simplifications -- If apply simplifications
    """
    lhs = sympy.Matrix(sympy.zeros(testfunc.shape[0],dofs.shape[0]) )
    for i in range(lhs.shape[0]):
        for j in range(lhs.shape[1]):
            lhs[i,j] = -sympy.diff(rhs[i,0], dofs[j,0])

            if do_simplifications:
                lhs[i,j] = sympy.simplify(lhs[i,j])

    return lhs

def Compute_RHS_and_LHS(functional, testfunc, dofs, do_simplifications=False):
    """ This computes the LHS matrix and the RHS vector

    Keyword arguments:
    functional -- The functional to derivate
    testfunc -- The test functions
    dofs -- The dofs vectors
    do_simplifications -- If apply simplifications
    """
    rhs = Compute_RHS(functional, testfunc, do_simplifications)
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications)
    return rhs,lhs

# Output functions
def _Indentation(indentation_level):
    "Returns the indentation string"
    return "    " * indentation_level

def _CodeGen(language, value):
    return  {
        "c"     : sympy.ccode,
        "python": sympy.pycode
    }[language](value)

def _VariableDeclaration(language, variable_name, variable_expression):
    """"Returns the variable declaration, without indentation nor suffix

    The expression must have been turned into code already
    """
    return  {
        "c"     : "const double {name} = {expr}",
        "python": "{name} = {expr}"
    }[language].format(name=variable_name, expr=variable_expression)

def _Suffix(language):
    "Returns the endline suffix"
    return  {
        "c"     : ";\n",
        "python": "\n"
    }[language]

def _ReplaceIndices(language, expression):
    """Replaces array access with underscored variable:
    For matrices: `variable[3,7]` becomes `variable_3_7`
    For vectors:  `variable[3]` becomes `variable_3`

    Depending on the language the accessors are chosen (`[]` vs. `()`)
    """
    #Matrices
    pattern = r"\[(\d+),(\d+)\]" if language == 'python' else r"\((\d+),(\d+)\)"
    replacement = r"_\1_\2"
    expression = re.sub(pattern, replacement, expression)

    # Vectors
    pattern = r"\[(\d+)\]" if language == 'python' else r"\((\d+)\)"
    replacement = r"_\1"
    expression = re.sub(pattern, replacement, expression)

    return expression


def OutputVector(r, name, language="python", initial_tabs=3, max_index=None, replace_indices=True, assignment_op="="):
    """ This method converts into text the RHS vector

    Keyword arguments:
    rhs -- The RHS vector
    name -- The name of the variables
    language -- The language of output
    initial_tabs -- The number of tabulations considered
    max_index -- DEPRECATED The maximum index
    replace_indices -- If the indixes must be replaced
    assignment_op -- The assignment operation
    """
    if max_index is not None:
        print("Warning: max_index parameter is deprecated in OutputVector")

    prefix = _Indentation(initial_tabs)
    suffix = _Suffix(language)
    fmt = prefix \
          + ("{var}[{i}]{op}{expr}" if language=="python" else "{var}({i}){op}{expr}") \
          + suffix

    outstring = str("")
    for i in range(r.shape[0]):
        expression = _CodeGen(language, r[i,0])
        outstring += fmt.format(var=name, i=i, op=assignment_op, expr=expression)

    if replace_indices:
        outstring = _ReplaceIndices(language, outstring)

    return outstring


def OutputMatrix(lhs, name, language, initial_tabs=3, max_index=None, replace_indices=True, assignment_op="="):
    """ This method converts into text the LHS matrix

    Keyword arguments:
    lhs -- The LHS matrix
    name -- The name of the variables
    language -- The language of output
    initial_tabs -- The number of tabulations considered
    max_index -- DEPRECATED The maximum index
    replace_indices -- If the indixes must be replaced
    assignment_op -- The assignment operation
    """
    if max_index is not None:
        print("Warning: max_index parameter is deprecated in OutputMatrix")

    prefix = _Indentation(initial_tabs)
    suffix = _Suffix(language)

    fmt = prefix \
          + ("{var}[{i},{j}]{eq}{expr}" if language == "python" else "{var}({i},{j}){eq}{expr}") \
          + suffix

    outstring = str("")
    for i in range(lhs.shape[0]):
        for j in range(lhs.shape[1]):
            expression = _CodeGen(language, lhs[i,j])
            outstring += fmt.format(var=name, i=i, j=j, op=assignment_op, expr=expression)

    if replace_indices:
        outstring = _ReplaceIndices(language, outstring)

    return outstring

def OutputSymbolicVariable(var, language="python", initial_tabs=None, max_index=None, replace_indices=True):
    """ This method converts into text the LHS matrix (only non-zero terms)
    Keyword arguments:
    var -- The variable to define symbolic
    language -- The language of output
    initial_tabs -- The number of tabulations considered
    max_index -- The maximum index
    replace_indices -- If the indixes must be replaced
    """
    if initial_tabs is not None:
        print("Warning: initial_tabs parameter is deprecated in OutputSymbolicVariable")
    if max_index is not None:
        print("Warning: max_index parameter is deprecated in OutputSymbolicVariable")

    outstring = _CodeGen(language, var) + _Suffix(language)

    if replace_indices:
        outstring = _ReplaceIndices(language, outstring)

    return outstring

def OutputSymbolicVariableAssignment(var, language, name, initial_tabs=3, max_index=None, replace_indices=True):
    """ This method converts into text the LHS matrix (only non-zero terms)

    Keyword arguments:
    var -- The variable to define symbolic
    language -- The language of output
    name -- The name of the variables
    initial_tabs -- The number of tabulations considered
    max_index -- DEPRECATED The maximum index
    replace_indices -- If the indixes must be replaced
    """

    if max_index is not None:
        print("Warning: max_index parameter is deprecated in OutputSymbolicVariable")

    prefix = _Indentation(initial_tabs)
    value = _CodeGen(language, var)
    expr = _VariableDeclaration(language, name, value)
    suffix = _Suffix(language)

    outstring = prefix + expr + suffix

    if replace_indices:
        _ReplaceIndices(language, outstring)

    return outstring

def _OutputX_CollectionFactors(A, name, language, initial_tabs, optimizations, replace_indices, assignment_op, output_func):
    """ This method collects the constants of the replacement for matrices, vectors, etc

    Keyword arguments:
    A -- The  factors
    name -- The name of the constant
    language -- The language of replacement
    initial_tabs -- The number of initial tabulations
    optimizations -- The level of optimizations
    replace_indices -- If the indixes must be replaced
    assignment_op -- The assignment operation
    output_func -- The output function. Must have the same signature as OutputMatrix and OutputVector
    """

    symbol_name = "c" + name
    A_factors, A_collected = sympy.cse(A, sympy.numbered_symbols(symbol_name), optimizations)
    A = A_collected[0] #overwrite lhs with the one with the collected components

    Acoefficient_str = str("")
    for factor in A_factors:
        varname = str(factor[0])
        value = factor[1]
        Acoefficient_str += OutputSymbolicVariableAssignment(value, language, varname, initial_tabs, None, replace_indices)

    A_out = Acoefficient_str + output_func(A, name, language, initial_tabs, None, replace_indices, assignment_op)
    return A_out


def OutputMatrix_CollectingFactors(A, name, language, initial_tabs=3, max_index=None, optimizations='basic', replace_indices=True, assignment_op="="):
    """ This method collects the constants of the replacement for matrices

    Keyword arguments:
    A -- The  factors
    name -- The name of the constant
    language -- The language of replacement
    initial_tabs -- The number of initial tabulations
    max_index -- DEPRECATED The max number of indexes
    optimizations -- The level of optimizations
    replace_indices -- If the indixes must be replaced
    assignment_op -- The assignment operation
    """
    if max_index is not None:
        print("Warning: max_index parameter is deprecated in OutputMatrix_CollectingFactors")

    return _OutputX_CollectionFactors(A, name, language, initial_tabs, optimizations, replace_indices, assignment_op, OutputMatrix)


def OutputVector_CollectingFactors(A, name, language, initial_tabs=3, max_index=None, optimizations='basic', replace_indices=True, assignment_op="="):
    """ This method collects the constants of the replacement for vectors

    Keyword arguments:
    A -- The  factors
    name -- The name of the constant
    language -- The language of replacement
    initial_tabs -- The number of initial tabulations
    max_index -- DEPRECATED The max number of indexes
    optimizations -- The level of optimizations
    replace_indices -- If the indixes must be replaced
    assignment_op -- The assignment operation
    """

    if max_index is not None:
        print("Warning: max_index parameter is deprecated in OutputVector_CollectingFactors")

    return _OutputX_CollectionFactors(A, name, language, initial_tabs, optimizations, replace_indices, assignment_op, OutputVector)
