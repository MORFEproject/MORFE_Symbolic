using MORFE_Symbolic
using SymPy
using Latexify

#include("MORFE_Symbolic.jl")

"""
Auxiliary function for latex printing. Creates non commutative monoms so that they stay in
the correct order when printing.
This function is primarely intended for use inside the module, and not for exportation.
"""
function non_commutative_monoms_for_latex_expressions(monom_number, monoms_exponents, vars)
    monom_string = ""
    for i = 1:length(vars)
        exp = monoms_exponents[monom_number][i]
        if exp != 0
            if exp == 1
                monom_string *= "z_$(i)*"
            else
                monom_string *= "z_$(i)^$(exp)*"
            end
        end
    end
    return sympy.Symbol(monom_string[1:end-1], commutative=False)
end

"""
Checks whether all of the terms in a sum have negative coefficients.
This function is primarely intended for use inside the module, and not for exportation.
"""
function all_negative_coeffs(expression)
    if expression.func == sympy.core.add.Add
        for arg in expression.args
            !sympy.core.function._coeff_isneg(arg) && return false
        end
    else
        !sympy.core.function._coeff_isneg(expression) && return false
    end
    return true
end

"""
Defines a consistent multiexponent ordering. For multiexponents m1 and m2,
if deg1 < deg2, then m1 < m2. If deg1 = deg2, then the comparison is made
looking at the degrees of each z_i in order. For example:
    - z1^3 * z2^4 * z3^1 > z1^3 * z2^3 * z3^3, because the z2 exponent is
    greater for the first monomial.
"""
function multiexponent_lt(x, y)
    sum(x) < sum(y) && return true
    sum(y) < sum(x) && return false

    for i in eachindex(x)
        x[i] > y[i] && return true    
        y[i] > x[i] && return false
    end
    return false
end

"""
Auxiliary function to transform a given polynomial expression in its latex code. Input arguments are:
    - vars: A symbolic vector containing the variables of the polynomial expression.
    - expr_coeff: The polynomial coefficients.
    - expr_monoms: A vector of tuples with each monom's exponents.
    - latex_output: The string to which the latex code is appended.
This function is primarely intended for use inside the module, and not for exportation.
"""
function latex_code_for_polynomial_expression(vars::Vector{Sym}, expr_coeffs::Vector{Sym}, expr_monoms, latex_output::String)
    ordering = sortperm(expr_monoms, lt=multiexponent_lt)
    expr_coeffs = expr_coeffs[ordering]
    expr_monoms = expr_monoms[ordering]

    for j in eachindex(expr_coeffs)
        sign_flag = false
        monom = non_commutative_monoms_for_latex_expressions(j, expr_monoms, vars)
        # This sign flag is used as a way to make the i stay outside of fraction numerators
        if all_negative_coeffs(expr_coeffs[j])
            expr_coeffs[j] = -expr_coeffs[j]
            sign_flag = true
        end
        if !(expr_coeffs[j]/im).has(Sym(im))
            expr_coeffs[j] = sympy.Mul(Sym(im), expr_coeffs[j]/im, evaluate=False)
        end
        # In the next lines the Sym("a") is necessary so that -i is not displayed as \left( -i \right)
        if expr_coeffs[j] == 1
            latex_result = latexify(sympy.Add(Sym("a"), monom, evaluate=False), cdot = false)[3:end-1]
        else
            latex_result = latexify(sympy.Add(Sym("a"), sympy.Mul(expr_coeffs[j], monom, evaluate=False), evaluate=False), cdot = false)[3:end-1]
        end
        if sign_flag
            latex_result = " -"*latex_result[3:end]
        end
        if j == 1 && !sign_flag
            latex_result = latex_result[3:end]
        end
        latex_output *= latex_result
    end
    return latex_output
end

"""
Function to output the reduced dynamics on latex format. Input arguments are:\\
    - DP: The parametrisation data structure.\\
    - aexp: A multiexponent_struct defining the numbering of the monoms in the parametrisation.\\
    - substitutions: An array of dictionaries (or of arrays of dictionaries) containing substitutions
                        to be made on the symbolic variables.\\
    - output_file: A file to output the results. If it is not passed, results are displayed
                    on the terminal.\\
    - file_mode: The opening mode of the file. If "a", the output is appended at the end.
                    If "w" overwrites existing files with the same name.
"""
function reduced_dynamics_latex_output(DP::parametrisation_struct, aexp::multiexponent_struct, substitutions, output_file = nothing, file_mode::String = "a")
    poly_from_expr = sympy.polys.polytools.poly_from_expr

    println("Printing reduced dynamics")

    zₜ = sympy.zeros(DP.n_rom,1)[:,1]
    z = [Sym("z_$(i)") for i=1:DP.n_rom]
    for i_ord=1:length(DP.f[1,:])
        monom = prod(z .^ aexp.mat[:,i_ord])
        for i_var=1:DP.n_rom
            substituted = mysub([DP.f[i_var,i_ord]],DP.subs[end:-1:1])
            for substitute in substitutions
                substituted = mysub(substituted, substitute)
            end
            zₜ[i_var:i_var] += substituted*monom
        end
    end
    
    latex_output = "\\begin{align}"
    for i in eachindex(zₜ)
        expr = poly_from_expr(zₜ[i], gens = z)
        expr_monoms = expr[1].monoms()
        expr_coeffs = expr[1].coeffs()
        latex_output *= "\n\\dot{z}_{$(i)} &="
        latex_output = latex_code_for_polynomial_expression(z, expr_coeffs, expr_monoms, latex_output)
        latex_output *= "\\\\"
    end
    latex_output = replace(latex_output, "I" => "\\mathit{i}")
    latex_output = latex_output[1:end-2] * "\n\\end{align}\n"
    
    if output_file === nothing
        println("Reduced dynamics:")
        println(latex_output)
    else
        open(output_file, file_mode) do file
            write(file, "Reduced dynamics:\n")
            write(file, latex_output)
        end 
    end
    println("Reduced dynamics printed")
end

"""
Function to output the nonlinear mappings on latex format. Input arguments are:\\
    - DP: The parametrisation data structure.\\
    - aexp: A multiexponent_struct defining the numbering of the monoms in the parametrisation.\\
    - substitutions: An array of dictionaries (or of arrays of dictionaries) containing substitutions
                        to be made on the symbolic variables.\\
    - output_file: A file to output the results. If it is not passed, results are displayed
                    on the terminal.\\
    - file_mode: The opening mode of the file. If "a", the output is appended at the end.
                    If "w" overwrites existing files with the same name.
"""
function nonlinear_mappings_latex_output(DP::parametrisation_struct, aexp::multiexponent_struct, substitutions, output_file = nothing, file_mode::String = "a")
    poly_from_expr = sympy.polys.polytools.poly_from_expr

    println("Printing nonlinear mappings")

    u = sympy.zeros(DP.n_aut,1)[:,1]
    z = [Sym("z_$(i)") for i=1:DP.n_rom]
    for i_ord=1:length(DP.W[1,:])
        monom = prod(z .^ aexp.mat[:,i_ord])
        for i_var=1:DP.n_aut
            substituted = mysub([DP.W[i_var,i_ord]],DP.subs[end:-1:1])
            for substitute in substitutions
                substituted = mysub(substituted, substitute)
            end
            u[i_var:i_var] += substituted*monom
        end
    end
    
    latex_output = "\\begin{align}"
    for i in eachindex(u)
        expr = poly_from_expr(u[i], gens = z)
        expr_monoms = expr[1].monoms()
        expr_coeffs = expr[1].coeffs()
        latex_output *= "\ny_{$(i)} &="
        latex_output = latex_code_for_polynomial_expression(z, expr_coeffs, expr_monoms, latex_output)
        latex_output *= "\\\\"
    end
    latex_output = replace(latex_output, "I" => "\\mathit{i}")
    latex_output = latex_output[1:end-2] * "\n\\end{align}\n"
    
    if output_file === nothing
        println("Nonlinear mappings:")
        println(latex_output)
    else
        open(output_file, file_mode) do file
            write(file, "Nonlinear mappings:\n")
            write(file, latex_output)
        end 
    end
    println("Nonlinear mappings printed")
end