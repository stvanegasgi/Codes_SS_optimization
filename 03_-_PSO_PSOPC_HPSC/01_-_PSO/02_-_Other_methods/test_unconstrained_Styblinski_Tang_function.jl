# Solved unconstrained Styblinski-Tang function in two variables with PSO flyback
# method
#
#                     n
# min  f(X) = (1/2) * Σ (xi^4 - 16xi^2 + 5xi)
#  X                 i=1
#
# X = [x, y]^T
#
# Bounds (or search domain):
# -5.0 <= xi <= 5.0
#
# Global minimun:
# X*      = [-2.903534, ..., -2.903534]^T
# -39.16617n <= f(X*) <= -39.16616n
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
# *(1) https://en.wikipedia.org/wiki/Test_functions_for_optimization
#
# -----------------------------------------------------------------------------

# ============================ packages =======================================

using Plots # for plots
ENV["GKSwstype"] = "100";

# ======================= PSO functions =======================================

# load the PSO flyback functions
include("./PSO_flyback_functions.jl");

# load the functions
include("../../../01_-_Optimization_models/unconstrained_Styblinski_Tang_function.jl");

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by f and f_c respectivily

num_var, _ = size(bounds); # number of variables

vbounds = [-2.0*ones(num_var)  2.0*ones(num_var)]; # velocity limits

opt_arg = nothing;  # optional arguments
num_particles = 50; # number of particles
k_max = 200;        # maximum iteration
ω   = [0.9, 0.6];   # inertia weight
c1  = [1.0, 1.0];   # acceleration constants
c2  = [1.0, 1.0];

values_f, X_g_value, population_data, Swarm_data, fun_evals, const_evals = PSO_flyback(f, f_c, bounds,
                                                                                       opt_arg,
                                                                                       vbounds,
                                                                                       num_particles,
                                                                                       k_max, ω, c1, c2);

# solution
solution    = Swarm_data.X_global;
f_solution  = Swarm_data.f_X_global;

println("\n\n=================================================================")
println("Unconstrained Styblinski Tang function --> PSO flyback")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(0:k_max, values_f, xlabel="Iterations k",
             ylabel="f(x)", label="")) # plot