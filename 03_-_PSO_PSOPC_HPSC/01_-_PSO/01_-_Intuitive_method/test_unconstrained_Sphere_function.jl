# Solved unconstrained Sphere function in multiple variables with PSO flyback method
#
#             n
# min  f(X) = Σ xi^2
#  X         i=1
#
# X = [x1, ..., xn]^T
#
# Bounds (or search domain):
# -∞ <= xi <= ∞
#
# Global minimun:
# X*      = [0.0, ..., 0.0]^T
# f(X*)   = 0
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
# *(1) https://en.wikipedia.org/wiki/Test_functions_for_optimization
#
# -----------------------------------------------------------------------------

# ======================= PSO functions =======================================

# load the PSO flyback functions
include("./PSO_flyback_functions.jl");

# load the functions
include("../../../01_-_Optimization_models/unconstrained_Sphere_function.jl");

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
println("Unconstrained Sphere function --> PSO flyback")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(0:k_max, values_f, xlabel="Iterations k",
             ylabel="f(x)", label="")) # plot
