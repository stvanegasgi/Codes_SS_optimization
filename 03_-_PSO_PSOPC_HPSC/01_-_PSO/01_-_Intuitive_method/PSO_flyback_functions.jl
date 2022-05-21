# Functions for solving optimization problem with the PSO flyback method.
#
# Optimization problem:
# min  f(X)     s.t  g_i(X) <= 0,   i = 1, 2, ..., m ---> constraints
#  X
#
# X = [x1, x2, ..., xn]^T
# xi^l <= xi <= xi^u,    i = 1, 2, ..., n ---> bounds
#
# =============================================================================
# DATE:    November 2021
# WHO:     Steven Vanegas Giraldo
# EMAIL:   stvanegasgi@unal.edu.co
# -----------------------------------------------------------------------------
# Universidad Nacional de Colombia - Sede Manizales
# =============================================================================
#
# References:
# *(1) Li LJ, Huang ZB, Liu F, Wu QH. A heuristic particle swarm optimizer for
#      optimization of pin connected structures.
#      Comput Struct 2007;85(7–8):340–9.
#
# *(2) He S, Prempain E, Wu QH. An improved particle swarm optimizer for
#      mechanical design optimization problems. Eng Optimiz 2004; 36(5):585–605.
#
# -----------------------------------------------------------------------------


# ============================= packages ======================================

using Distributions; # for probabilistic distributions

# ============================= data structures ===============================

# data struct for each particle
mutable struct Particle
    X       ::Vector{Float64}   # actual position
    X_prev  ::Vector{Float64}   # previous position
    f_X     ::Float64           # function value of actual position
    V       ::Vector{Float64}   # actual velocity
    X_best  ::Vector{Float64}   # best position
    f_X_best::Float64           # function value of best position
end

# data struct of swarm
mutable struct Swarm
    number_particles::Int64             # number of particles in the swarm
    best_particle   ::Int64             # best particle in swarm
    X_global        ::Vector{Float64}   # best position in swarm
    f_X_global      ::Float64           # function value of X global in swarm
    V_max_comp      ::Vector{Float64}   # maximum value in velocity components
    V_min_comp      ::Vector{Float64}   # minimum value in velocity components
end

# ============================= auxiliar functions ============================

"""
Function to determine if X lies between the boundaries

X_feasible = X_feasible_bounds(X_test, X_min, X_max)

Parameters:

    X_test  (Array):    variable to test
                        -- array type -- Vector{Float64}

    X_min   (Array):    lower bounds of variables
                        -- array type -- Vector{Float64}

    X_max   (Array):    upper bounds of variables
                        -- array type -- Vector{Float64}

Returns:

    X_feasible      (Bool):     flag indicates that X lies between boundaries
                                (true/yes) or (false/not)
"""
function X_feasible_bounds(X_test::Vector{Float64},
                           X_min ::Vector{Float64},
                           X_max ::Vector{Float64})

#   check bounds
    check_upper_bound = all(X_test .<= X_max);
    check_lower_bound = all(X_test .>= X_min);

#   X is between the bounds?
    X_feasible = (check_upper_bound && check_lower_bound);

    return X_feasible
end


"""
Function to determine if X complies with constraints

    X_feasible = X_feasible_constraints(X_test, constraint_fun, optional_arg)

Parameters:

    X_test             (Array):     variable to test
                                    -- array type -- Vector{Float64}

    constraint_fun  (Function):     function that evaluate all the constraints
                                    of the problem in the form g(X) <= 0
                                    -------------------------------------------
                                    constraint_fun(X, opt_arg)
                                    X:       array with the value to compute
                                    opt_arg: optional argument

    optional_arg         (Any):     optional argument for the function of
                                    constraints

Returns:

    X_feasible      (Bool):     flag indicates that X complies all the
                                constraints (true/yes) or (false/no)
"""
function X_feasible_constraints(X_test        ::Vector{Float64},
                                constraint_fun::Function,
                                optional_arg  ::Any)

#   counter constraint evaluations
    global counter_constraint;
    counter_constraint += 1;

#   all the constraints are minus o equal 0, gi(X) <= 0 (constraints vector)
    X_feasible = all(constraint_fun(X_test, optional_arg) .<= 0);

    return X_feasible
end


"""
Computes a linear variation of parameters

    value = linear_variation(k, k_max, initial_value, final_value)

Parameters:

    k               (Int):      actual iteration

    k_max           (Int):      maximum iteration for the method

    initial_value (Float):      initial value

    final_value   (Float):      final value

Returns:

    value         (Float):      value for the iteration k
"""
function linear_variation(k             ::Int64,
                          k_max         ::Int64,
                          initial_value ::Float64,
                          final_value   ::Float64)

#   linear function of inertia weight
    value = ((final_value - initial_value)/(k_max - 1))*(k - 1) + initial_value;

    return value
end


"""
Inits the variables of the particle

    particle = init_particle(X_min, X_max, constraint_fun, fun,
                             opt_arg, Vmin, Vmax)

Parameters:

    X_min              (Array):   array with the minimum bounds for each
                                  component of X
                                  -- array type -- Vector{Float64}

    X_max              (Array):   array with the maximum bounds for each
                                  component of X
                                  -- array type -- Vector{Float64}

    constraint_fun  (Function):   function that evaluate all the constraints
                                  of the problem in the form g(X) <= 0
                                  -------------------------------------------
                                  constraint_fun(X, opt_arg)
                                  X:       array with the value to compute
                                  opt_arg: optional argument

    fun             (Function):   function that computes the value of the
                                  objective function to optimize, with the input
                                  parameters as follows:
                                  ---------------------------------------------
                                  fun(X, opt_arg)
                                  X:       array with the value to compute
                                  opt_arg: optional argument

    optional_arg         (Any):   optional argument for the function of
                                  constraints

    Vmin               (Array):   minimum bounds for the velocity
                                  -- array type -- Vector{Float64}

    Vmax               (Array):   maximum bounds for the velocity
                                  -- array type -- Vector{Float64}

Returns:

    particle          (Struct):   struct particle (propieties)
"""
function init_particle(X_min            ::Vector{Float64},
                       X_max            ::Vector{Float64},
                       constraint_fun   ::Function,
                       fun              ::Function,
                       opt_arg          ::Any,
                       Vmin             ::Vector{Float64},
                       Vmax             ::Vector{Float64})

#   defines the counters variables
    global counter_funtion
    global counter_constraint

    check_constraint = false;   # X is not feasible in constraints

    local X_part                # local scope out the while loop
    while check_constraint == false

#       X particle between the bounds (for each component)
        X_part = rand.(Uniform.(X_min, X_max));

#       counter constraint function evaluations
        counter_constraint += 1;

#       is feasible X in constraint?
        check_constraint = X_feasible_constraints(X_part,
                                                  constraint_fun,
                                                  opt_arg);
    end

#   clip the velocities for each particle
    V_part = rand.(Uniform.(Vmin, Vmax));

#   the particle data is initialize
    particle = Particle(deepcopy(X_part), # actual position (X_part)
                        deepcopy(X_part), # previous position
                        Inf,              # function value of X_part
                        deepcopy(V_part), # actual velocity
                        deepcopy(X_part), # best position
                        Inf          );   # function value of best position

#   counter function evaluations
    counter_funtion += 1;

#   value function in actual position a best position of particle
    particle.f_X      = fun(particle.X, opt_arg);
    particle.f_X_best = deepcopy(particle.f_X);

    return particle
end


"""
Update the particle velocity, between the bounds

    particle = update_velocity(dim_X, particle, Swarm_i, ω, c1, c2, Vmin, Vmax)

Parameters:

    dim_X           (Int): dimension of X

    particle     (Struct): struct particle (propieties)

    Swarm_i      (Struct): struct of swarm (propieties)

    ω             (Float): inertia weight in the actual iteration

    c1            (Float): acceleration constant in the actual iteration,
                           for the cognitive velocity

    c2            (Float): acceleration constant in the actual iteration,
                           for the social velocity

    Vmin          (Array): minimum bounds for the velocity
                           -- array type -- Vector{Float64}

    Vmax          (Array): maximum bounds for the velocity
                           -- array type -- Vector{Float64}

Returns:

    particle          (Struct):   struct particle (propieties)
"""
function update_velocity(dim_X      ::Int64,
                         particle   ::Particle,
                         Swarm_i    ::Swarm,
                         ω          ::Float64,
                         c1         ::Float64,
                         c2         ::Float64,
                         Vmin       ::Vector{Float64},
                         Vmax       ::Vector{Float64})

#   uniform random sequence U(0, 1)
    r1, r2 = rand(Uniform(), dim_X), rand(Uniform(), dim_X);

    V = deepcopy(particle.V); # actual particle velocity

#   updating the particles velocity (new velocity)
    particle.V = ω*V +                                    # inertia
                 c1*r1.*(particle.X_best  - particle.X) + # cognitive
                 c2*r2.*(Swarm_i.X_global - particle.X);  # social

#   clip under the limits
    particle.V = clamp.(particle.V, Vmin, Vmax);

    return particle
end


# ============================= main function =================================

"""
Function that applies the PSO flyback method. Reference [1]

fun_g_value, X_g_value, population, Swarm_prop = PSO_flyback(fun, c_fun, 
                                                             Xbounds, opt_arg,
                                                             Vbounds,
                                                             num_particles,
                                                             k_max, ω, c1, c2)

Parameters:

    fun            (Function):  function that computes the value of the
                                objective function to optimize, with the input
                                parameters as follows:
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

    Xbounds           (Array):  array that contain all the lower and upper
                                bounds for each component of vector X, for
                                n components
                                -----------------------------------------------
                                Xbounds = [x1^l      x1^u;
                                           x2^l      x2^u;
                                           ...       ... ;
                                           xn^l      xn^u]
                                -- array type -- Array{Float64, 2}

    opt_arg             (Any):  optional argument for the objective function or
                                constraints function. The structure depends on
                                how the functions fun and c_fun were built

    Vbounds           (Array):  array that contain the lower and upper bounds
                                for the components of the velocity vector
                                -----------------------------------------------
                                Vbounds = [Vmin_1  Vmax_1;
                                           Vmin_2  Vmax_2;
                                             ...    ...  ;
                                           Vmin_n  Vmax_n]
                                for default:
                                Vmin   = -abs.(Xmax - Xmin);
                                Vmax   =  abs.(Xmax - Xmin);
                                -- array type -- Array{Float64, 2}

    num_particles       (Int):  number of particles. For default 20

    k_max               (Int):  maximum number of iterations. For default 200

    ω                 (Array):  array with the initial and final values of the
                                inertia weight. The value of ω changes linearly
                                with the iterations.
                                For default ω = [0.8, 0.8] --> constant
                                -----------------------------------------------
                                ω = [ω_initial, ω_final]
                                -- array type -- Array{Float64, 1}

    c1                (Array):  array with the initial and final values of the
                                constant c1. The value of c1 changes linearly
                                with the iterations. Acceleration constant
                                for cognitive velocity.
                                For default c1 = [2.0, 2.0] --> constant
                                -----------------------------------------------
                                c1 = [c1_initial, c1_final]
                                -- array type -- Array{Float64, 1}

    c2                (Array):  array with the initial and final values of the
                                constant c2. The value of c2 changes linearly
                                with the iterations. Acceleration constant
                                for social velocity.
                                For default c2 = [2.0, 2.0] --> constant
                                -----------------------------------------------
                                c2 = [c2_initial, c2_final]
                                -- array type -- Array{Float64, 1}

Returns:

    fun_g_value       (Array):  array with all value of objetive function of
                                the best global position in the swarm for each
                                iteration
                                -- array type -- Array{Float64, 1}

    X_g_value         (Array):  array with all the best global solution of the
                                swarm for each iteration. The column k is a X
                                global for the iteration k
                                -- array type -- Array{Float64, 2}

    population        (Array):  array with the propieties of all particles,
                                each element is a data struct for each particle

    Swarm_prop       (Struct):  struct with the propieties of the swarm

    counter_funtion     (Int):  counter objective function evaluations

    counter_constraint  (Int):  counter constraint evaluations
"""
function PSO_flyback(fun            ::Function,
                     c_fun          ::Function,
                     Xbounds        ::Array{Float64,2},
                     opt_arg        ::Any=nothing,
                     Vbounds        ::Array{Float64,2}=nothing,
                     num_particles  ::Int64=20,
                     k_max          ::Int64=200,
                     ω              ::Vector{Float64}=[0.8, 0.8],
                     c1             ::Vector{Float64}=[2.0, 2.0],
                     c2             ::Vector{Float64}=[2.0, 2.0])

# ========================= lecture variables =================================

#   variables for lower and upper limits
    LOWER = 1;
    UPPER = 2;

# ========================= some variables ====================================

#   counter for objective function evaluation and restrictions function
    global counter_funtion = 0;
    global counter_constraint = 0;

    Xmin  = Xbounds[:, LOWER];         # minimun bounds for X
    Xmax  = Xbounds[:, UPPER];         # maximum bounds for X
    dim_X = length(Xbounds[:, LOWER]); # dimensions of X

    if !=(Vbounds, nothing)
        Vmin = Vbounds[:, LOWER]; # minimun bounds for velocity
        Vmax = Vbounds[:, UPPER]; # maximum bounds for velocity
    else # default values
        Vmin = -abs.(Xmax - Xmin);
        Vmax =  abs.(Xmax - Xmin);
    end

#   minimum and maximum parameters of method
    ω_initial  =  ω[LOWER];    ω_final  =  ω[UPPER];
    c1_initial = c1[LOWER];    c1_final = c1[UPPER];
    c2_initial = c2[LOWER];    c2_final = c2[UPPER];

#   initialize some propieties of swarms
    Swarm_prop = Swarm(num_particles,   # number of particles
                       1,               # number of best particle
                       zeros(dim_X),    # best position in swarm
                       Inf,             # function value of best position in swarm
                       deepcopy(Vmax),  # maximum value in velocity components
                       deepcopy(Vmin)); # minimun value in velocity components

#   array with all struct of particles is initialize
    population = Array{Particle, 1}(undef, num_particles);

# ============== initialize the positions and velocities of each particle =====

    k = 1; # initialize iteration counter
    for i_particle in 1:num_particles

#       init the particle variables, between the bounds
        particle = init_particle(Xmin, Xmax, c_fun, fun, opt_arg, Vmin, Vmax);

#       is a best global solution?
        if particle.f_X < Swarm_prop.f_X_global
            Swarm_prop.best_particle = i_particle;
            Swarm_prop.f_X_global    = deepcopy(particle.f_X);
            Swarm_prop.X_global      = deepcopy(particle.X);
        end

        population[i_particle] = deepcopy(particle); # new particle
    end

# ============================ main loop ======================================

    fun_g_value = zeros(k_max+1);         # historical function value (global)
    X_g_value   = zeros(dim_X, k_max+1);  # historical X (global) for each
#                                           iteration

#   save the best value function global and X solution for each iteration (k=0)
    fun_g_value[k]  = Swarm_prop.f_X_global;
    X_g_value[:, k] = deepcopy(Swarm_prop.X_global);

#   print the initial solution (k = 0)
    println("\n\nK [$(0)]: f(X_global) = $(Swarm_prop.f_X_global) ")
    println("X_global = $(Swarm_prop.X_global)")

#   main loop
    while k <= k_max

#       proposal X_global and it value function for the iteration
        X_global_i     = deepcopy(Swarm_prop.X_global);
        fun_X_global_i = deepcopy(Swarm_prop.f_X_global);
        particle_X_g_i = deepcopy(Swarm_prop.best_particle);

#       for each particle
        for i_particle in 1:num_particles

            particle        = population[i_particle]; # actual particle
            particle.X_prev = deepcopy(particle.X);   # update the previous X

#           computes the parameter value for iteration k
            ω_value  = linear_variation(k, k_max, ω_initial,  ω_final);
            c1_value = linear_variation(k, k_max, c1_initial, c1_final);
            c2_value = linear_variation(k, k_max, c2_initial, c2_final);

#           update the particle velocity
            particle = update_velocity(dim_X, particle, Swarm_prop,
                                       ω_value, c1_value, c2_value, Vmin, Vmax);

            particle.X += particle.V; # updating the new position

#           is feasible in bounds and constraints?
            check_bounds      = X_feasible_bounds(particle.X, Xmin, Xmax);
            check_constraints = X_feasible_constraints(particle.X, c_fun,
                                                       opt_arg);

#           is not feasible in bounds and constraints?
            if ~(check_bounds && check_constraints)
#               return X_(k) to X_(k-1), flyback mechanism
                particle.X = deepcopy(particle.X_prev);
            end

            counter_funtion += 1; # conunter function evaluations
#           value function in actual position
            particle.f_X = fun(particle.X, opt_arg);

#           is a best solution for the particle?
            if particle.f_X < particle.f_X_best
                particle.X_best   = deepcopy(particle.X);
                particle.f_X_best = particle.f_X;
            end

#           is a best global solution? (it's a proposal for next iteration)
            if particle.f_X < fun_X_global_i
                particle_X_g_i = i_particle;
                fun_X_global_i = particle.f_X;
                X_global_i     = deepcopy(particle.X);
            end
        end

#       update the global values
        Swarm_prop.best_particle = particle_X_g_i;
        Swarm_prop.f_X_global    = fun_X_global_i;
        Swarm_prop.X_global      = deepcopy(X_global_i);

#       print the results each 5 iterations
        if (k % 5) == 0 || k == 1
            println("\nK [$(k)]: f(X_global) = $(Swarm_prop.f_X_global) ")
            println("X_global = $(Swarm_prop.X_global)")
        end

#       save the best value function global and X solution
        fun_g_value[k+1]  = Swarm_prop.f_X_global;
        X_g_value[:, k+1] = deepcopy(Swarm_prop.X_global);

        k += 1 # next iteration
    end

#   print the final results
    println("\n\nFinal results:")
    println("K [$(k_max)]: f(X_global) = $(Swarm_prop.f_X_global) "*
            "X_global = $(Swarm_prop.X_global)")

    return fun_g_value, X_g_value, population, Swarm_prop, counter_funtion, counter_constraint
end
