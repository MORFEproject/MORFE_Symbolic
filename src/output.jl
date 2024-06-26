using SymPy
using Latexify

"""
Auxiliary function for latex printing. Creates non commutative monoms so that they stay in
the correct order when printing.
This function is primarely intended for use inside the module, and not for exportation.
"""
function non_commutative_monoms_for_latex_expressions(monom_number, monoms_exponents, normal_coordinate)
    if sum(monoms_exponents[monom_number]) == 0
        return 1
    end
    monom_string = ""
    for i = 1:length(monoms_exponents[1])
        exp = monoms_exponents[monom_number][i]
        if exp != 0
            if exp == 1
                monom_string *= normal_coordinate * "_$(i)*"
            else
                monom_string *= normal_coordinate * "_$(i)^$(exp)*"
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
This function is primarely intended for use inside the module, and not for exportation.
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
function latex_code_for_polynomial_expression(expr_coeffs::Vector{Sym}, expr_monoms, latex_output::String, normal_coordinate)
    ordering = sortperm(expr_monoms, lt=multiexponent_lt)
    expr_coeffs = expr_coeffs[ordering]
    expr_monoms = expr_monoms[ordering]

    for j in eachindex(expr_coeffs)
        sign_flag = false
        monom = non_commutative_monoms_for_latex_expressions(j, expr_monoms, normal_coordinate)
        # This sign flag is used as a way to make the i stay outside of fraction numerators
        if all_negative_coeffs(expr_coeffs[j])
            expr_coeffs[j] = -expr_coeffs[j]
            sign_flag = true
        end
        if !(simplify(expr_coeffs[j]/im)).has(Sym(im))
            expr_coeffs[j] = sympy.Mul(Sym(im), simplify(expr_coeffs[j]/im), evaluate=False)
        end
        # In the next lines the Sym("a") is necessary so that -i is not displayed as \left( -i \right)
        if monom == 1
            latex_result = latexify(sympy.Add(Sym("a"), expr_coeffs[j], evaluate=False), cdot = false, safescripts = true)
        elseif expr_coeffs[j] == 1
            latex_result = latexify(sympy.Add(Sym("a"), monom, evaluate=False), cdot = false, safescripts = true)
        else
            latex_result = latexify(sympy.Add(Sym("a"), sympy.Mul(expr_coeffs[j], monom, evaluate=False), evaluate=False), cdot = false, safescripts = true)
        end
        if latex_result[2] == 'a'
            latex_result = latex_result[3:end-1]
        else
            latex_result = " + " * latex_result[2:end-5]
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
function reduced_dynamics_latex_output(DP::parametrisation_struct, aexp::multiexponent_struct, output_file = nothing;
                                       file_mode::String = "a", normal_coordinate = 'z', result = "complex")
    poly_from_expr = sympy.polys.polytools.poly_from_expr

    println("Printing reduced dynamics")

    zₜ = sympy.zeros(DP.n_rom,1)[:,1]
    z = [Sym("$(normal_coordinate)_$(i)") for i=1:DP.n_rom]
    for i_ord=1:length(DP.f[1,:])
        monom = prod(z .^ aexp.mat[:,i_ord])
        if result == "complex"
            zₜ += DP.f[:,i_ord]*monom
        elseif result == "real"
            zₜ += DP.fr[:,i_ord]*monom
        else 
            throw(ArgumentError("Result type should be either complex or real!"))
        end 
    end
    
    latex_output = "\\begin{align}"
    for i in eachindex(zₜ)
        expr = poly_from_expr(zₜ[i], gens = z)
        expr_monoms = expr[1].monoms()
        expr_coeffs = expr[1].coeffs()
        latex_output *= "\n\\dot{$(normal_coordinate)}_{$(i)} &="
        latex_output = latex_code_for_polynomial_expression(expr_coeffs, expr_monoms, latex_output, normal_coordinate)
        latex_output *= "\\\\"
    end
    latex_output = replace(latex_output, "I" => "\\mathit{i}")
    latex_output = latex_output[1:end-2] * "\n\\end{align}\n"
    
    if output_file === nothing
        if result == "complex"
            println("Reduced dynamics:")
        else
            println("Realified reduced dynamics:")
        end 
        println(latex_output)
    else
        open(output_file, file_mode) do file
            if result == "complex"
                write(file, "Reduced dynamics:\n")
            else
                write(file, "Realified reduced dynamics:\n")
            end
            write(file, latex_output)
        end 
    end
    println("Reduced dynamics printed\n")
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
function nonlinear_mappings_latex_output(DP::parametrisation_struct, aexp::multiexponent_struct, output_file = nothing;
                                         file_mode::String = "a", normal_coordinate = 'z', result = "complex")
    poly_from_expr = sympy.polys.polytools.poly_from_expr

    println("Printing nonlinear mappings")

    u = sympy.zeros(DP.n_full,1)[:,1]
    z = [Sym("$(normal_coordinate)_$(i)") for i=1:DP.n_rom]
    for i_ord=1:length(DP.W[1,:])
        monom = prod(z .^ aexp.mat[:,i_ord])
        if result == "complex"
            u += DP.W[:,i_ord]*monom
        elseif result == "real"
            u += DP.Wr[:,i_ord]*monom
        elseif result == "modal"
            u += DP.Wmodal[:,i_ord]*monom
        else 
            throw(ArgumentError("Result type should be either complex, real or modal!"))
        end 
    end
    
    latex_output = "\\begin{align}"
    for i in eachindex(u)
        expr = poly_from_expr(u[i], gens = z)
        expr_monoms = expr[1].monoms()
        expr_coeffs = expr[1].coeffs()
        if result == "modal"
            latex_output *= "\ny_{$(i)} &="
        else
            if i <= DP.n_osc
                latex_output *= "\nu_{$(i)} &="
            elseif i <= 2*DP.n_osc
                latex_output *= "\nv_{$(i-DP.n_osc)} &="
            else
                latex_output *= "\nr_{$(i-2*DP.n_osc)} &="
            end 
        end
        latex_output = latex_code_for_polynomial_expression(expr_coeffs, expr_monoms, latex_output, normal_coordinate)
        latex_output *= "\\\\"
    end
    latex_output = replace(latex_output, "I" => "\\mathit{i}")
    latex_output = latex_output[1:end-2] * "\n\\end{align}\n"
    
    if output_file === nothing
        if result == "complex"
            println("Nonlinear mappings:")
        elseif result == "real"
            println("Realified nonlinear mappings:")
        else
            println("Nonlinear mappings to modal coordinates:")
        end 
        println(latex_output)
    else
        open(output_file, file_mode) do file
            if result == "complex"
                write(file, "Nonlinear mappings:\n")
            elseif result == "real"
                write(file, "Realified nonlinear mappings:\n")
            else
                write(file, "Nonlinear mappings to modal coordinates:\n")
            end 
            write(file, latex_output)
        end 
    end
    println("Nonlinear mappings printed\n")
end

"""
Function to output the backbone obtained by polar realification. Only to be used when the problem
is parametrised by a single master mode. Input arguments are:\\
    - omega_rho: An array containing the monomials defining the backbone when it is seen as a polynomial in rho;\\
    - output_file: A file to output the results. If it is not passed, results are displayed
                    on the terminal.\\
    - file_mode: The opening mode of the file. If "a", the output is appended at the end.
                    If "w" overwrites existing files with the same name.
"""
function backbone_output(omega_rho, output_file = nothing; file_mode::String = "a")
    poly_from_expr = sympy.polys.polytools.poly_from_expr

    println("Printing backbone")
    ρ = symbols("ρ", real=true)
    latex_output = "\\begin{align}"

    sum = 0
    for i in eachindex(omega_rho)
        sum += omega_rho[i]
    end
    expr = poly_from_expr(sum, gens = ρ)
    expr_monoms = expr[1].monoms()
    expr_coeffs = expr[1].coeffs()
    latex_output *= "\n\\omega &="
    latex_output = latex_code_for_polynomial_expression(expr_coeffs, expr_monoms, latex_output, "ρ")
    latex_output *= "\\\\"
    latex_output = replace(latex_output, "I" => "\\mathit{i}")
    latex_output = latex_output[1:end-2] * "\n\\end{align}\n"

    if output_file === nothing
        println("Backbone:")
        println(latex_output)
    else
        open(output_file, file_mode) do file
            write(file, "Backbone:\n")
            write(file, latex_output)
        end 
    end
    println("Backbone printed\n")
end

"""
Function to output the nonlinear damping ratio obtained by polar realification. Only to be used when the problem
is parametrised by a single master mode. Input arguments are:\\
    - xi_rho: An array containing the monomials defining the damping when it is seen as a polynomial in rho;\\
    - output_file: A file to output the results. If it is not passed, results are displayed
                    on the terminal.\\
    - file_mode: The opening mode of the file. If "a", the output is appended at the end.
                    If "w" overwrites existing files with the same name.
"""
function nonlinear_damping_output(xi_rho, output_file = nothing; file_mode::String = "a")
    poly_from_expr = sympy.polys.polytools.poly_from_expr

    println("Printing backbone")
    ρ = symbols("ρ", real=true)
    latex_output = "\\begin{align}"

    sum = 0
    for i in eachindex(xi_rho)
        sum += xi_rho[i]
    end
    expr = poly_from_expr(sum, gens = ρ)
    expr_monoms = expr[1].monoms()
    expr_coeffs = expr[1].coeffs()
    latex_output *= "\n\\xi_{NL} &="
    latex_output = latex_code_for_polynomial_expression(expr_coeffs, expr_monoms, latex_output, "ρ")
    latex_output *= "\\\\"
    latex_output = replace(latex_output, "I" => "\\mathit{i}")
    latex_output = latex_output[1:end-2] * "\n\\end{align}\n"

    if output_file === nothing
        println("Nonlinear damping:")
        println(latex_output)
    else
        open(output_file, file_mode) do file
            write(file, "Nonlinear damping:\n")
            write(file, latex_output)
        end 
    end
    println("Nonlinear damping printed\n")
end

"""
Function to output the physical amplitude obtained by polar realification. Only to be used when the problem
is parametrised by a single master mode. Input arguments are:\\
    - ampli_rho: An array containing the monomials defining the amplitude when it is seen as a polynomial in rho;\\
    - output_file: A file to output the results. If it is not passed, results are displayed
                    on the terminal.\\
    - file_mode: The opening mode of the file. If "a", the output is appended at the end.
                    If "w" overwrites existing files with the same name.
"""
function physical_amplitudes_output(ampli_rho, output_file = nothing; file_mode::String = "a")
    poly_from_expr = sympy.polys.polytools.poly_from_expr

    println("Printing physical amplitudes")
    ρ = symbols("ρ", positive=true)
    latex_output = "\\begin{align}"

    sum = 0
    for i in eachindex(ampli_rho)
        sum += ampli_rho[i]
    end
    expr = poly_from_expr(sum, gens = ρ)
    expr_monoms = expr[1].monoms()
    expr_coeffs = expr[1].coeffs()
    latex_output *= "\nu_{max} &="
    latex_output = latex_code_for_polynomial_expression(expr_coeffs, expr_monoms, latex_output, "ρ")
    latex_output *= "\\\\"
    latex_output = replace(latex_output, "I" => "\\mathit{i}")
    latex_output = latex_output[1:end-2] * "\n\\end{align}\n"

    if output_file === nothing
        println("Physical amplitudes:")
        println(latex_output)
    else
        open(output_file, file_mode) do file
            write(file, "Physical amplitudes:\n")
            write(file, latex_output)
        end 
    end
    println("Physical amplitudes printed\n")
end

"""
Function to write .nb files in order to make a link with Mathematica. Always prints a file with all of the
relevant variables and may optionally also print separate files with the reduced dynamics and nonlinear 
mappings. Input arguments are:\\
    - DP: The parametrisation data structure.\\
    - aexp: A multiexponent_struct defining the numbering of the monoms in the parametrisation.\\
    - directory: The directory at which the files are going to be saved.\\
    - file_basename: The base name for the files. It will be appended with either "_variables", 
    "_reduced_dynamics" or "_nonlinear_mappings".\\
    - print_reduced_dynamics: logic flag to determine if the reduced dynamics should be printed.\\
    - print_nonlinear_mappings: logic flag to determine if the nonlinear mappings should be printed.
"""
function Mathematica_output(DP::parametrisation_struct, aexp::multiexponent_struct, directory, file_basename;
    print_reduced_dynamics = true, print_nonlinear_mappings = true)
    mathematica_code = sympy.printing.mathematica.mathematica_code

    println("Mathematica output started")

    isdir(directory) || mkpath(directory)
    path = joinpath(directory, file_basename * "_variables.nb")
    open(path, "w") do file
        write(file, "naut = $(DP.n_aut);\n")
        write(file, "nfull = $(DP.n_full);\n")
        write(file, "nsets = $(DP.n_sets);\n")
        write(file, "order = $(DP.order);\n")
        for i in eachindex(DP.subs)
            for key in keys(DP.subs[i])
                if key != 0
                    math_code = mathematica_code(key) * " = FullSimplify[" * mathematica_code(DP.subs[i][key]) * "];\n"
                    math_code = replace(math_code, "_" => "")
                    write(file, math_code) 
                end
            end
        end
        write(file, "monoms = ConstantArray[0,$(DP.n_sets)];\n")
        for i_set = 1:DP.n_sets
            z = [Sym("z$(i)") for i=1:DP.n_rom]
            monom = prod(z .^ aexp.mat[:,i_set])
            monom = mathematica_code(monom)
            write(file, "monoms[[$(i_set)]] = " * monom * ";\n")
        end
        if print_reduced_dynamics
            write(file, "f = ConstantArray[0,{$(DP.n_rom),$(DP.n_sets)}];\n")
            for i_var = 1:DP.n_rom
                for i_set = 1:DP.n_sets
                    math_code = "f[[$(i_var),$(i_set)]] = FullSimplify[" * mathematica_code(DP.f[i_var,i_set]) * "];\n"
                    math_code = replace(math_code, "_" => "")
                    write(file, math_code)
                end
            end
        end
        if print_nonlinear_mappings
            write(file, "W = ConstantArray[0,{$(DP.n_full),$(DP.n_sets)}];\n")
            for i_var = 1:DP.n_full
                for i_set = 1:DP.n_sets
                    math_code = "W[[$(i_var),$(i_set)]] = FullSimplify[" * mathematica_code(DP.W[i_var,i_set]) * "];\n"
                    math_code = replace(math_code, "_" => "")
                    write(file, math_code)
                end
            end
        end
    end

    if print_reduced_dynamics
        path = joinpath(directory, file_basename * "_reduced_dynamics.nb")
        open(path, "w") do file
            write(file, "fres = ConstantArray[0, naut];\n")
            write(file, "For[ivar = 1, ivar <= naut, ivar++,\n")
            write(file, "   For[iset = 1, iset <= nsets, iset++,\n")
            write(file, "       fres[[ivar]] += f[[ivar, iset]]*monoms[[iset]];\n")
            write(file, "   ];\n")
            write(file, "];\n")
        end
    end

    if print_nonlinear_mappings
        path = joinpath(directory, file_basename * "_nonlinear_mappings.nb")
        open(path, "w") do file
            write(file, "Wres = ConstantArray[0, nfull];\n")
            write(file, "For[ivar = 1, ivar <= nfull, ivar++,\n")
            write(file, "   For[iset = 1, iset <= nsets, iset++,\n")
            write(file, "       Wres[[ivar]] += W[[ivar, iset]]*monoms[[iset]];\n")
            write(file, "   ];\n")
            write(file, "];\n")
        end
    end
    println("Mathematica output finished")
end