# Functions for solving optimization problem with the Subset Simulation method
#
# Optimization problem:
# max  f(X)     s.t  g_i(X) <= 0,   i = 1, 2, ..., m ---> constraints
#  X
#
# X = [x1, x2, ..., xn]^T
# xi^l <= xi <= xi^u,    i = 1, 2, ..., n ---> bounds
#
# =============================================================================
# DATE:    March 2022
# WHO:     Steven Vanegas Giraldo
# EMAIL:   stvanegasgi@unal.edu.co
# -----------------------------------------------------------------------------
# Universidad Nacional de Colombia - Sede Manizales
# =============================================================================
#
# References:
# *(1) Li, H.-S., Au, S.-K (2010). Desing optimization using subset simulation
#      algorithm. Structural Safety, 32(6), 384-392.
#
# -----------------------------------------------------------------------------

# ============================= packages ======================================

using Distributions; # for probabilistic distributions
using Statistics;    # basic statistics functionality

# ============================= auxiliar functions ============================

"""
Function to compute constraint violation function vi, equation 6, Ref. [1]

    v_i = const_violation_function(value_gi)

Parameters:

    value_gi        (Float):    value of the constraint gi(x)

Returns:

    v_i             (Float):    value of the constraint violation function
                                equation (6) in [1]
                                -----------------------------------------------
                                         ( 0        if  gi(x) <= 0
                                v_i(x) = |
                                         ( gi(x)    if  gi(x) >  0
"""
function const_violation_function(value_gi)

    # Equation 6, Ref. [1]
    if value_gi > 0     # for gi(xs) > 0  (xs don't complies with the restrictions)
        return value_gi
    else                # for gi(xs) <= 0 (xs complies with the restrictions)
        return 0
    end
end

"""
Function to compute the constraint fitness function for a sample Xs

    Fcon_Xs     =   F_con(Xs                   ::Vector{Float64},
                          constraint_function  ::Function,
                          opt_arg              ::Any=nothing)

Parameters:

    Xs                  (Array):    sample to test Xs = [x1, x2, ..., xN]^T

    constraint_function (Function): function that computes the value of the all
                                    constraints of optimization problem, the
                                    constraints must have the form gi(X) <= 0, the
                                    array must contain all the values gi(X) <= 0.
                                    With the input parameters as follows:
                                    -------------------------------------------
                                    constraint_function(X, opt_arg)
                                    X:       array with the value to compute
                                    opt_arg: optional argument

    opt_arg             (Any):      optional argument for the objective function or
                                    constraints function. The structure depends on
                                    how the function constraint_function were built

Returns:

    Fcon                (Float):    value of Fcon(xs) = -max(v_i), equation (7)
                                    in [1]
"""
function F_con(Xs                    ::Vector{Float64},
               constraint_function   ::Function,
               opt_arg               ::Any=nothing)

#   evaluate the constrains vector [g1(x), g2(x), ..., gL(x)]^T
    constrains = constraint_function(Xs, opt_arg);

#   evaluate the constraint violation function, equation (6) in [1]
    v_i = const_violation_function.(constrains);

#   constraint fitness function equation (7) in [1]
    return -max(v_i...)
end

"""
Function to compute constraint fitness function for a sample Xs

    fun_x_sort, Fc_x_sort, samples_sort = double_criterion_sort(fun_x  ::Vector{Float64},
                                                                Fc_x   ::Vector{Float64},
                                                                samples::Array{Float64,2})

Parameters:

    fun            (Function):  function that computes the value of the
                                objective function to optimize, with the input
                                parameters, as follows:
                                -----------------------------------------------
                                fun(X, opt_arg)
                                X:       array with the value to compute
                                opt_arg: optional argument

    Fc_x              (Array):  array with the values of constraint fitness
                                function of the samples (unsorted)
                                Fc_x = [Fcon(x1), Fcon(x2), ..., Fcon(xN)]^T

    samples           (Array):  array with the samples, each row is a dimension
                                and each column is a sample, N samples and n
                                components, matrix (n x N)
                                -----------------------------------------------
                                samples =  [x1^1   x1^2        x1^N;
                                            x2^1   x2^2        x2^N;
                                            ...    ...         ... ;
                                            xn^1   xn^2        xn^N]

Returns:

    fun_x_sort        (Array):  array with the values of objective function sorted
                                by the double_criterion

    Fc_x_sort         (Array):  array with the values of constraint fitness function
                                sorted by the double_criterion

    samples_sort      (Array):  array with the samples sorted by the double criterion
                                the first column are the best values in terms of
                                Fcon() and value of objective function (complies
                                for samples with Fcon(-) = 0)
"""
function double_criterion_sort(fun_x    ::Vector{Float64},
                               Fc_x     ::Vector{Float64},
                               samples  ::Array{Float64,2})

#   first sort the constrains fitness function, Fcon_x, equation (8) [1]
    ind_1     = sortperm(Fc_x, rev=true); # find the index
    Fc_x_sort = deepcopy(Fc_x[ind_1]);    # first the best sample

#   verify the samples that satisfy Fcon = 0
    check_Fcon = findall(values -> values.== 0, Fc_x_sort);

#   final index sort
    final_index = deepcopy(ind_1);

#   only sort with objective function the feasible solutions
    if length(check_Fcon) != 0 # there exist feasibles samples
        fun_x_sort = deepcopy(fun_x[ind_1]); # sort the objective functions values

#       second sort by objective function only complies Fcon = 0 (feasible samples)
#       best the max of objective function value
        ind_2 = sortperm(fun_x_sort[check_Fcon], rev=true)

#       sort the feasible samples with the objective function
        final_index[check_Fcon] = final_index[check_Fcon][ind_2];
    end

#   ultimate sort (the double-criterion) first the best sample
    fun_x_sort   = deepcopy(fun_x[final_index]);      # objetive function
    Fc_x_sort    = deepcopy(Fc_x[final_index]);       # Fcon
    samples_sort = deepcopy(samples[:, final_index]); # samples

    return fun_x_sort, Fc_x_sort, samples_sort
end



# ============================= main function =================================

"""
Subset simulation (ss) algorithm for design optimization. Reference [1]

    x_optimal, h_optimal, f_samples_k_level, hk_k_level, Fconk_k_level = ss_optimization(fun     ::Function,
                                                                                         c_fun   ::Function,
                                                                                         N       ::Int64,
                                                                                         bounds  ::Array{Float64,2},
                                                                                         opt_arg ::Any=nothing,
                                                                                         ε       ::Float64)

Parameters:

    fun            (Function):  function that computes the value of the
                                objective function to optimize, with the input
                                parameters, as follows:
                                -----------------------------------------------
                                fun(X, opt_arg)
                                X:       array with the value to compute
                                opt_arg: optional argument

    c_fun          (Function):  function that computes the value of the all
                                constraints of optimization problem, the
                                constraints must have the form gi(X) <= 0, the
                                array must be contain all the values gi(X) <= 0.
                                With the input parameters as follows:
                                -----------------------------------------------
                                c_fun(X, opt_arg)
                                X:       array with the value to compute
                                opt_arg: optional argument

    N                   (Int):  number of samples for the method ss

    bounds            (Array):  array that contains the lower and upper
                                bounds for each component of vector X, for
                                n components
                                -----------------------------------------------
                                bounds =  [x1^l      x1^u;
                                           x2^l      x2^u;
                                           ...       ... ;
                                           xn^l      xn^u]

    opt_arg             (Any):  optional argument for the objective function or
                                constraints function. The structure depends on
                                how the functions fun and c_fun were built

    ε                 (Float):  convergence criterion, equation (13) in [1]
                                max(abs.(σ_k./(bounds[:, upper] - bounds[:, lower])))

Returns:

    x_optimal       (Array):    array with all the best x for each level
                                -----------------------------------------------
                                x_optimal =  [x1_opt_1      x1_opt_end;
                                              x2_opt_1      x2_opt_end;
                                               ...             ...    ;
                                              xn_opt_1      xn_opt_end] 

    h_optimal       (Array):    array with all value of objetive function of
                                the best x for each level
                                -----------------------------------------------
                                h_optimal = [h_1, h_2, ..., h_end]^T

    samples_k_level (Array):    array with the samples for each k-level. Each
                                element is an array with all the samples. The
                                rows in the array are components and each column
                                is a sample

    f_samples_k_level (Array):  array with the objective function value for each
                                set of samples per k-level is a function of
                                samples_k_level. Each element is an array with
                                the function values for the samples in the k-level

    hk_k_level      (Array):    array with the funtion threshold values for
                                each k-level

    Fconk_k_level   (Array):    array with the Fcon threshold values for each
                                k-level
"""
function ss_optimization(fun            ::Function,
                         c_fun          ::Function,
                         N              ::Int64,
                         bounds         ::Array{Float64,2},
                         opt_arg        ::Any=nothing,
                         ε              ::Float64=1e-6)

# ========================= lecture variables =================================

#   indices for upper and lower limit
    LOWER = 1;
    UPPER = 2;

# ========================= some variables ====================================

#   means of the boundaries (for each row) center of domain
    μ_i = (bounds[:, LOWER] + bounds[:, UPPER])/2;

#   three sigma limits, equation (12) [1] (3 standard deviation from the mean)
    σ_i = abs.(bounds[:, UPPER] - bounds[:, LOWER])/6;

#   dimension of the variables X = [x1, x2, ..., xn]^T
    dimension = length(μ_i);   

# ========================= compute P(F1) =====================================

#   prepares samples matrix (row=component, column=sample)
    XS = zeros(dimension, N);

#   generate the sample of f(x; μ, σ, xl, xu) for each component
    for component in 1:dimension

        xl_i, xu_i  = bounds[component, :]; # lower and upper limit

#       truncated normal distribution (artificial PDF), equation (11) [1]
        trunc_normal = Truncated(Normal(μ_i[component], σ_i[component]), xl_i, xu_i);

#       uniform distribution (N samples for component) U ~ [0,1] for samples
        U = rand(Uniform(0,1), N);

#       generate the samples from the artificial PDF for the component (cdf^-1(U))
        xs_i = quantile(trunc_normal, U);

#       include the sample
        XS[component, :] = xs_i; 
    end

    Fcon_x = zeros(N); # prepared constraint fitness function for each sample
    h_x    = zeros(N); # evaluated objective function for each sample

#   evaluate the samples xs_i (objective function h(x) and constraint fitness
#   function F_con(x) equation (7) [1]) for each sample
    for sample = 1:N
#       objetive function for the sample
        h_x[sample]    = fun(XS[:, sample], opt_arg);

#       constrain fitness function for the sample
        Fcon_x[sample] = F_con(XS[:, sample], c_fun, opt_arg);
    end

#   sort the samples for double-criterion in [1]
    h_x, Fcon_x, XS = double_criterion_sort(h_x, Fcon_x, XS);

#   initial stage
    k = 1;

    σ_k = std(XS, dims=2); # standard desviation for the samples for components

#   convergence criterion, equation (13) in [1]
    σ_k_bounds = abs.(σ_k./(bounds[:, UPPER] - bounds[:, LOWER]));

#   p_1 = 0.5 (initial level probability) ϵ [0, 1]
    p_k = 0.5;

#   save x_optimal and h(x_optimal) for the k-level
    x_optimal = [deepcopy(XS[:, 1])];
    h_optimal = [h_x[1]];

#   prepared the variable for samples and objective function per k level
    samples_k_level   = [deepcopy(XS)];  # is a list of arrays
    f_samples_k_level = [deepcopy(h_x)]; # is a array

#   variables to save the thresholds per k-level 
    hk_k_level    = Array{Float64, 1}();
    Fconk_k_level = Array{Float64, 1}();

    while ε < maximum(σ_k) # convergence criterion

#       print the results for each level
        println("x =  ", x_optimal[:, end])
        println("h(x) =  ", h_optimal[end])
        println("level =  ", k, "\n")

#       initial calculations for the k-th simulation level
        Nc = ceil(Int64, p_k * N); # number of chains
        Ns = ceil(Int64, 1/p_k);   # number of samples per chain

        if N != Nc * Ns # Nc and Ns are positive integers?
            error("Nc or Ns are not a positive numbers. Stop simulation\n");
        end

        Ns = Ns - 1; # fit the samples

#       prepared the seeds for the k simulation level
        seeds      = deepcopy(XS[:, 1:Nc]);  # the seeds
        h_seeds    = deepcopy(h_x[1:Nc]);    # objective function of the seeds
        Fcon_seeds = deepcopy(Fcon_x[1:Nc]); # constraints fitness function of the seeds
        h_k  = h_x[Nc];    # the threshold h_k for the k-level
        Fc_k = Fcon_x[Nc]; # the threshold Fc_k for the k-level

#       save the thresholds
        push!(hk_k_level, h_k);
        push!(Fconk_k_level, Fc_k);

        for chain in 1:Nc # number of Markov chain

            seed_i       = deepcopy(seeds[:, chain]); # seed with n component
            h_seeds_i    = h_seeds[chain];    # objetive function of the seed
            Fcon_seeds_i = Fcon_seeds[chain]; # constraint fitness function of the seed

            for sample_chain in 1:Ns # for Ns samples per chain

#               a candidate X'
                X_p = zeros(dimension); # prepared a sample X' from p*(-|Xk)

                for component_i in 1:dimension # for each component of the sample
#                                                                               OJO ES MAS FACIL GENERARLO DE LA PDF TRUNCADA con Truncated
                    xl_i, xu_i = bounds[component_i, :];

#                   the proposal PDF p*(-|Xk), N(Xk, σ_sample)
                    proposal_pdf = Normal(seed_i[component_i], σ_k[component_i]);

#                   in some case the sample could be outer the boundaries
                    local x_component;
                    while true # break if the sample is between the boundaries
#                       generate sample for the component i for the proposal PDF
                        x_component = rand(proposal_pdf);

#                       the sample could be outer the boundaries
                        if xl_i <= x_component <= xu_i
                            break;
                        end
                    end

#                   compute the probability of acceptation r
#                   truncated normal distribution (artificial PDF), equation (11) [1]
                    trunc_normal = Truncated(Normal(μ_i[component_i], σ_i[component_i]), xl_i, xu_i);

#                   relation f_j(X')/f_j(Xk)
                    r = pdf(trunc_normal, x_component)/pdf(trunc_normal, seed_i[component_i]);

                    r = min(1, r); # acceptance probability

#                   acceptation the component i of sample or is Xk
                    if rand(Uniform(0, 1)) < r  # X_k+1 = X'
                        X_p[component_i] = x_component;         # accept
                    else                        # X_k+1 = Xk
                        X_p[component_i] = seed_i[component_i]; # reject
                    end
                end

#               evaluate a objetive function and the constraint fitness function
#               for the candidate
                h_candidate    = fun(X_p, opt_arg);
                Fcon_candidate = F_con(X_p, c_fun, opt_arg);

                idx = Nc + (chain - 1)*Ns + sample_chain; # index

#               determinate the candidate X' is in the event k-level with Fk
                if Fc_k == 0 # satisfy the constraint {h >= hk}

                    if h_candidate >= h_k # accept the candidate in the chain
#                       save the candidate, objective function and his Fcon
                        XS[:, idx]  = deepcopy(X_p);
                        h_x[idx]    = h_candidate;
                        Fcon_x[idx] = Fcon_candidate;
                    else # reject the candidate and set the seed
#                       save the seed, objective function and his Fcon
                        XS[:, idx]  = deepcopy(seed_i);
                        h_x[idx]    = h_seeds_i;
                        Fcon_x[idx] = Fcon_seeds_i;
                    end
                else # Fck < 0 {Fcon >= Fck}

                    if Fcon_candidate >= Fc_k # accept the candidate in the chain
#                       save the candidate, objective function and his Fcon
                        XS[:, idx]  = deepcopy(X_p);
                        h_x[idx]    = h_candidate;
                        Fcon_x[idx] = Fcon_candidate;
                    else # reject the candidate and set the seed
#                       save the seed, objective function and his Fcon
                        XS[:, idx]  = deepcopy(seed_i);
                        h_x[idx]    = h_seeds_i;
                        Fcon_x[idx] = Fcon_seeds_i;
                    end
                end

#               define the news seed and objective function of the seed for the chain
                seed_i       = deepcopy(XS[:, idx]);
                h_seeds_i    = h_x[idx];
                Fcon_seeds_i = Fcon_x[idx];
            end
        end

#       sort the samples for double-criterion in [1]
        h_x, Fcon_x, XS = double_criterion_sort(h_x, Fcon_x, XS);

        σ_k = std(XS, dims=2); # standard desviation for the new samples for components

#       convergence criterion, equation (13) in [1]
        σ_k_bounds = abs.(σ_k./(bounds[:, UPPER] - bounds[:, LOWER]));

#       defines the level probability for the next level, section 4.3 [1]
        σ_i_hat = maximum(σ_k); # criterion for the level probability

#       criterion from section 4.3 [1]
        if σ_i_hat < 0.1
            p_k = 0.2;
            if σ_i_hat < 0.01
                p_k = 0.1
            end
        end

#       save x_optimal and h(x_optimal) for the k-level
        x_optimal = [x_optimal deepcopy(XS[:, 1])];
        h_optimal = [h_optimal; h_x[1]];

#       save the samples and their objective function per k-level
        push!(samples_k_level, deepcopy(XS));
        push!(f_samples_k_level, deepcopy(h_x));

        k = 1 + k; # nest k_level
    end

    return x_optimal, h_optimal, samples_k_level, f_samples_k_level, hk_k_level, Fconk_k_level
end
