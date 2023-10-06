function create_multiexp(n::Int,o::Int)
    #= 
    a is a dictionary containing both forward and inverse mappings 
    from a certain multiexponent set [a1;a2;a3..;ap] 
    to its "position" i and order p
    but also its opposite: from [p i] to the set
    moreover, it also contains the number of sets for a certain p

    get_a is a function that interrogates the dictionary
    if given a vector [a1;a2;..;ap] will give the 1x2 matrix [p i]
    if given a 1x2 matrix [p i] will give the vector [a1;a2;..;ap]
    if given an integer p will give the number of vectors [a1;a2;..;ap] of that order 
    =#
    mltexp = collect(multiexponents(n,0));
    MltExp = Matrix{Int}(undef, n,1);
    for in = 1:n;        MltExp[in,1] = mltexp[1][in];    end
    a = Dict(mltexp[1] => 1)
    a = merge(a, Dict([0 1] =>mltexp[1]))
    a = merge(a,Dict(0 =>1))
    n_sets = 1
    for p = 1:o
        mltexp = collect(multiexponents(n,p));
        mltexp_mat = Matrix{Int}(undef, n,length(mltexp))
        for in = 1:n;  for im = 1:length(mltexp);    mltexp_mat[in,im] = mltexp[im][in];    end;end
        MltExp = [MltExp mltexp_mat];
        a = merge(a,Dict(MltExp[:,n_sets+i] => n_sets+i for i = 1:length(mltexp)))
        a = merge(a,Dict([p i] => MltExp[:,n_sets+i] for i = 1:length(mltexp)))
        a = merge(a,Dict(p =>length(mltexp)))
        n_sets = n_sets+length(mltexp);
    end
    function get_a(a_set)
        return get(a,a_set,"-")
    end
    return get_a,a,MltExp,n_sets
end

function partitions_two(p::Int)
    # partitions gives a list of all the 2 integers whose sum is p
    # except the trivial one [p,0] and [0,p]
    #
    #
    # the part vector will be for example:
    #= 
    
    partitions_two(5)    = 
    4-element Vector{Vector{Int64}}:
    [4, 1]
    [3, 2]
    [2, 3]
    [1, 4]
    
    =#
    part = collect(partitions(p,2))
    if iseven(p)
        for k = length(part)-1:-1:1
            part = [part;[part[k][end:-1:1]]]
        end
    else
        for k = length(part):-1:1
            part = [part;[part[k][end:-1:1]]]
        end
    end
    return part
end

"""
Struct to manage multiexponents. It is constituted of the following:
    - dic: A dictionary containing pairs of (key, value). There are three kinds of keys:
        - [deg local_index]      : Used to access the monomial of degree deg with number local_index.
                                    For example [4, 3] accesses the 3rd monomial of degree 4.
                                    The returned value is an array, in the format of next item in this list.
        - [deg_1, ..., deg_n_rom]: Used to access of the monomial z_1^deg_1 * ... * z_n^deg_n_rom,
                                    where n_rom represents the number of master coordinates.
                                    The returned value is an integer, the global index of the monomial.
        - deg                    : Returns an integer, the number of multiexponents of degree deg.
    - get: A function to get values from dic.
    - mat: An n_rom x n_sets matrix, in which each column is a multiexponent.
    - n_sets: The number of monoms with n_rom variables up to order o.
"""
mutable struct multiexponent_struct
    get
    dic::Dict{Any,Any}
    mat::Matrix{Int}
    n_sets::Int
end

"""
Initializes the multiexponent struct for a model with n_rom master coordinates and
parametrization up to order o.
"""
function init_multiexponent_struct(n_rom::Int,o::Int)
    g,d,v,n_sets = create_multiexp(n_rom,o)
    return multiexponent_struct(g,d,v,n_sets)
end