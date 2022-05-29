# Functions for solving optimization problem with the crow search algorithm (CSA)
#
#
# Optimization problem:
# min  f(X)     s.t  g_i(X) <= 0,   i = 1, 2, ..., m ---> constraints
#  X
#
# X = [x1, x2, ..., xn]^T
# xi^l <= xi <= xi^u,    i = 1, 2, ..., n ---> bounds
#
# =============================================================================
# DATE:    May 2022
# WHO:     Steven Vanegas Giraldo
# EMAIL:   stvanegasgi@unal.edu.co
# -----------------------------------------------------------------------------
# Universidad Nacional de Colombia - Sede Manizales
# =============================================================================
#
# References:
# * Askarzadeh, A. (2016). A novel metaheuristic method for solving constrained
#   engineering optimization problems: crow search algorithm.
#   Computers & Structures, 169, 1-12.
#
# -----------------------------------------------------------------------------

# ============================= packages ======================================

using Distributions; # for probabilistic distributions

# ============================= data structures ===============================

# data struct for each crow
mutable struct Crow
    x   ::Vector{Float64} # actual position in actual iteration
    f_x ::Float64         # function value of actual position
    m   ::Vector{Float64} # hiding place memorazed in actual iteration
    f_m ::Float64         # function value of hiding place memorazed
end;

# ============================ auxiliar functions =============================

"""
Function to determine if x complies with boundaries

    x_feasible = x_feasible_boundaries(x_test, x_min, x_max)

Parameters:

    x_test  (Array):     variable to test
                         -- array type -- Vector{Float64}

    x_min   (Array):    array with the minimum bounds of the variables
                        -------------------------------------------------------
                        -- array type -- Vector{Float64}
            
    x_max   (Array):    array with the maximus bounds of the variables
                        -------------------------------------------------------
                        -- array type -- Vector{Float64}

Returns:

    x_feasible      (Bool):     flag indicates that x complies all the
                                boundaries (true/yes) or (false/no)
"""
function x_feasible_boundaries(x_test::Vector{Float64},
                               x_min ::Vector{Float64},
                               x_max ::Vector{Float64})

#   check bounds
    check_upper_bound = all(x_test .<= x_max);
    check_lower_bound = all(x_test .>= x_min);

#   x_test is between the bounds?
    x_feasible = (check_upper_bound && check_lower_bound);

    return x_feasible # is feasible in constraints?
end;


"""
Function to determine if x complies with constraints

    x_feasible = x_feasible_constraints(x_test, constraint_fun, optional_arg)

Parameters:

    x_test             (Array):     variable to test
                                    -- array type -- Vector{Float64}

    constraint_fun  (Function):     function that evaluate all the constraints
                                    of the problem in the form g(x) <= 0
                                    -------------------------------------------
                                    constraint_fun(x, opt_arg)
                                    x:       array with the value to compute
                                    opt_arg: optional argument

    optional_arg         (Any):     optional argument for the function of
                                    constraints

Returns:

    x_feasible      (Bool):     flag indicates that x complies all the
                                constraints (true/yes) or (false/no)
"""
function x_feasible_constraints(x_test        ::Vector{Float64},
                                constraint_fun::Function,
                                optional_arg  ::Any)

    global num_eval_constraints; # global variable
    num_eval_constraints += 1; # constraints evaluations

#   all the constraints are minus or equal 0, gi(x) <= 0 (constraints vector)
    x_feasible = all(constraint_fun(x_test, optional_arg) .<= 0);

    return x_feasible # is feasible in constraints?
end;


"""
Function to inits the crow's data

    crow_struct = init_crow(f, f_c, X_min, X_max, opt_arg)

Parameters:

    f    (Function):    function that evaluate the objective function of the
                        problem
                        -------------------------------------------------------
                        f(x, opt_arg)
                        x:       array with the value to compute
                        opt_arg: optional argument

    f_c  (Function):    function that evaluate all the constraints of the
                        problem in the form g(x) <= 0
                        -------------------------------------------------------
                        f_c(x, opt_arg)
                        x:       array with the value to compute
                        opt_arg: optional argument

    X_min   (Array):    array with the minimum bounds of the variables
                        -------------------------------------------------------
                        -- array type -- Vector{Float64}

    X_max   (Array):    array with the maximus bounds of the variables
                        -------------------------------------------------------
                        -- array type -- Vector{Float64}

    opt_arg   (Any):    optional argument for the function of constraints and
                        objective function

Returns:

    crow_struct (Struct):   flag indicates that x complies all the
                            constraints (true/yes) or (false/no)
"""
function init_crow(f      ::Function,
                   f_c    ::Function,
                   X_min  ::Vector{Float64},
                   X_max  ::Vector{Float64},
                   opt_arg::Any)

#   counter by global variables
    global num_eval_funtion;

    is_feasible = false; # is a feasible position?

    local crow_x; # local scope out the while loop
    while is_feasible == false

        crow_x = rand.(Uniform.(X_min, X_max)); # a random position

#       is a feasible position?
        is_feasible = x_feasible_constraints(crow_x, f_c, opt_arg);
    end

    num_eval_funtion += 1; # objective function evaluations
    f_crow_x = f(crow_x, opt_arg); # eval the crow positions

#   crow's data
    crow = Crow(deepcopy(crow_x),   # initial position
                f_crow_x,           # objective function of position
                deepcopy(crow_x),   # hiding place
                f_crow_x);          # objective function of hiding place

    return crow
end;


# ============================ main function ==================================

"""
Function to aplied the crow seach algorithm (CSA)

    x_opti_iter, fun_x_opti_iter, num_eval_funtion, num_eval_constraints = CSA(fun, fun_c, x_min, x_max, N, iter_max, fl, AP, optional_arg)

Parameters:

    fun     (Function):  function that evaluate the objective function of the
                         problem
                         ------------------------------------------------------
                         fun(x, optional_arg)
                         x:            array with the value to compute
                         optional_arg: optional argument

    fun_c   (Function):  function that evaluate all the constraints of the
                         problem in the form g(x) <= 0
                         ------------------------------------------------------
                         fun_c(x, optional_arg)
                         x:            array with the value to compute
                         optional_arg: optional argument

    x_min      (Array):  array with the minimum bounds for each component of x
                         -- array type -- Vector{Float64}

    x_max      (Array):  array with the maximum bounds for each component of x
                        -- array type -- Vector{Float64}

    N          (Int64):  number of crows in the flock. For default 100 crow

    iter_max   (Int64):  maximus iterations in the algorithm. For default 300
                         iterations

    fl       (Float64):  flight length of crows. For default 2

    AP       (Float64):  awareness probability of crows. For default 0.1

    optional_arg (Any):  optional argument for the function of constraints and
                         objective function

Returns:

    x_opti_iter          (Array): array with all the best global solution of the
                                  flock for each iteration. The column k is a x
                                  global for the iteration k
                                  -- array type -- Array{Float64, 2}

    fun_x_opti_iter      (Array): array with all value of objetive function of
                                  the best global position in the flock for each
                                  iteration
                                  -- array type -- Array{Float64, 1}

    num_eval_funtion     (Int64): number of objective function evaluations

    num_eval_constraints (Int64): number of constraints evaluations
"""
function CSA(fun            ::Function,
             fun_c          ::Function,
             x_min          ::Vector{Float64},
             x_max          ::Vector{Float64},
             N              ::Int64=100,
             iter_max       ::Int64=300,
             fl             ::Float64=2.0,
             AP             ::Float64=0.1,
             optional_arg   ::Any=nothing)

#   dimension of the problem
    dim = length(x_min);

#   variables of the best solutions in the flock and his function value for each iteration
    x_opti_iter     = zeros(dim, iter_max + 1);
    fun_x_opti_iter = zeros(iter_max + 1);

#   temporal variable for the best solution in the flock for each iteration
    x_opti   = zeros(dim);
    f_x_opti = Inf64;

#   variables for evaluations objective function and constraints
    global num_eval_funtion     = 0;
    global num_eval_constraints = 0;

#   array with all struct of crow is initialize
    flock = Array{Crow, 1}(undef, N);

#   for each crow initialize positions and memories
    for crow_i = 1:N
#       init the feasible crow position in the boundaries
        crow_struct_i = init_crow(fun, fun_c, x_min, x_max, optional_arg);

        if crow_struct_i.f_m < f_x_opti # best solution in the flock is searched

            f_x_opti = crow_struct_i.f_m;
            x_opti   = deepcopy(crow_struct_i.m);
        end

        flock[crow_i] = deepcopy(crow_struct_i); # save in the flock
    end

    iter = 1; # begin the iterations (initial values)

#   save best solution in the flock
    x_opti_iter[:, iter]  = deepcopy(x_opti);
    fun_x_opti_iter[iter] = f_x_opti;

#   print the initial solution (k = 0)
    println("\n\nK [$(0)]: f(X*) = $(fun_x_opti_iter[iter]) ")
    println("X* = $(x_opti_iter[:, iter])")

#   main loop
    while iter <= iter_max

        for i = 1:N # for each crow

            Crow_i = flock[i]; # the actual crow

            r_j = rand(Uniform()); # number r_j ~ U(0,1)

            if r_j >= AP # (awareness probability)

#               randomly choose one of the crows to follow
                random_crow = rand(DiscreteUniform(1, N));
                Crow_j      = flock[random_crow];

                r_i = rand(Uniform()); # number r_i ~ U(0,1)

#               compute a new position
                x_new = Crow_i.x + r_i * fl * (Crow_j.m - Crow_i.x);
            else
#               a random position in the seach space
                x_new = rand.(Uniform.(x_min, x_max));
            end

#           is a feasible position?
            check_bounds          = x_feasible_boundaries(x_new, x_min, x_max);
            check_constraints     = x_feasible_constraints(x_new, fun_c, optional_arg);

            if check_bounds && check_constraints # only update the feasible positions

#               update the new feasible position and objective function
                Crow_i.x   = deepcopy(x_new);
                num_eval_funtion += 1;
                Crow_i.f_x = fun(Crow_i.x, optional_arg);

#               the hiding place is renovated
                if Crow_i.f_x < Crow_i.f_m

                    Crow_i.m   = deepcopy(x_new); # new hiding place
                    Crow_i.f_m = Crow_i.f_x;      # new objective function of
#                                                   new hiding place
                end
                flock[i] = deepcopy(Crow_i); # save the crow's struct
            end

#           best solution in the flock?
            if Crow_i.f_m < f_x_opti
                f_x_opti = Crow_i.f_m;
                x_opti   = deepcopy(Crow_i.m);
            end
        end

        iter += 1; # new iteration

#       save best solution in the flock
        x_opti_iter[:, iter]  = deepcopy(x_opti);
        fun_x_opti_iter[iter] = f_x_opti;

#       print the results each 5 iterations
        if (iter % 5) == 0
            println("\n\nK [$(iter)]: f(X*) = $(fun_x_opti_iter[iter]) ")
            println("X* = $(x_opti_iter[:, iter])")
        end
    end

    return x_opti_iter, fun_x_opti_iter, num_eval_funtion, num_eval_constraints
end;
