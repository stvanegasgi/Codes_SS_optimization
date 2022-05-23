# Solve Gomez and Levy function (modified) with the PSO flyback method
#
#
# min  f(X) = 4x^2 - 2.1x^4 + (1/3)x^6 + xy - 4y^2 + 4y^4
#  X
#
# X = [x, y]^T
#
# s.t: -sin(4 pi x) + 2 sin(2 pi y)^2 - 1.5 <= 0
#
# Bounds:
# -1.0 <= x <= 0.75
# -1.0 <= y <= 1.0
#
# Global minimun:
# X*     = [0.08984201, -0.7126564]^T
# f(X*)  = -1.031628453
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
# * https://en.wikipedia.org/wiki/Test_functions_for_optimization
#
# -----------------------------------------------------------------------------

# ============================ packages =======================================

using Plots # for plots
ENV["GKSwstype"] = "100";

# ======================= PSO functions =======================================

# load the PSO flyback functions
include("./PSO_flyback_functions.jl");

# load the functions
include("../../../01_-_Optimization_models/constrained_Gomez_Levy_function.jl");

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by f and f_c respectivily

num_var, _ = size(bounds); # number of variables

vbounds = [-2.0*ones(num_var)  2.0*ones(num_var)]; # velocity limits

opt_arg = nothing;  # optional arguments
num_particles = 30; # number of particles
k_max = 150;        # maximum iteration
ω   = [0.7, 0.5];   # inertia weight
c1  = [2.0, 2.0];   # acceleration constants
c2  = [2.0, 2.0];

values_f, X_g_value, population_data, Swarm_data, fun_evals, const_evals = PSO_flyback(f, f_c, bounds,
                                                                                       opt_arg,
                                                                                       vbounds,
                                                                                       num_particles,
                                                                                       k_max, ω, c1, c2);

# solution
solution    = Swarm_data.X_global;
f_solution  = Swarm_data.f_X_global;

println("\n\n=================================================================")
println("Constrained Gomez Levy function --> PSO flyback")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(0:k_max, values_f, xlabel="Iterations k",
             ylabel="f(x)", label="")) # plot
