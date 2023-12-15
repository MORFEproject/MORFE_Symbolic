function create_gen_vec(str,n)
    z = sympy.zeros(n,1)[:,1]
    for i = 1:n
        z[i] = symbols(str*string(i))
    end
    return z
end

function create_pos_vec(str,n)
    z = sympy.zeros(n,1)[:,1]
    for i = 1:n
        z[i] = symbols(str*string(i),positive = true)
    end
    return z
end

function mysimp(input)
    exp = reshape(input,prod(size(input)))
    for i = 1:length(exp)
        exp[i] = simplify(exp[i])
    end
    return reshape(exp,size(input))
end

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