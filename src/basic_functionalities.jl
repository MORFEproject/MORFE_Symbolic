"""
Creates a generic symbolic vector where each entry is given by str indexed by an integer, up to size n.
"""
function create_gen_vec(str,n)
    z = sympy.zeros(n,1)[:,1]
    for i = 1:n
        z[i] = symbols(str*"_"*string(i))
    end
    return z
end

"""
Creates a real symbolic vector where each entry is given by str indexed by an integer, up to size n.
"""
function create_real_vec(str,n)
    z = sympy.zeros(n,1)[:,1]
    for i = 1:n
        z[i] = symbols(str*"_"*string(i),real = true)
    end
    return z
end

"""
Creates a positive symbolic vector where each entry is given by str indexed by an integer, up to size n.
"""
function create_pos_vec(str,n)
    z = sympy.zeros(n,1)[:,1]
    for i = 1:n
        z[i] = symbols(str*"_"*string(i),positive = true)
    end
    return z
end

"""
Creates a real matrix where each entry is given by str indexed by an integer, up to size n.
"""
function create_real_matrix(str,n)
    T = sympy.zeros(n,n)
    for i = 1:n
        for j = 1:n
            T[i,j] = symbols(str*"_"*string(i)*string(j),real = true)
        end
    end
    return T
end

"""
Creates a positive matrix where each entry is given by str indexed by an integer, up to size n.
"""
function create_pos_matrix(str,n)
    T = sympy.zeros(n,n)
    for i = 1:n
        for j = 1:n
            T[i,j] = symbols(str*"_"*string(i)*string(j),positive = true)
        end
    end
    return T
end

"""
Simplifies the input
"""
function mysimp(input)
    exp = reshape(input,prod(size(input)))
    for i = 1:length(exp)
        exp[i] = simplify(exp[i])
    end
    return reshape(exp,size(input))
end

"""
Performs substitutions on the input. subs should be a vector of dictionaries, each one having
the item to substitute as key and the substitution as value. 
"""
function mysub(input,sub)
    exp = reshape(input,prod(size(input)))
    for i = 1:length(exp)
        for j = 1:length(sub)
            exp[i] = exp[i] |> subs(sub[j])
        end
        for j = 1:length(sub)
            exp[i] = expand(exp[i]) |> subs(sub[j])
            exp[i] = simplify(exp[i])
        end
    end
    return reshape(exp,size(input))
end