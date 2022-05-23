# Solve Rosenbrock function with the PSO flyback method (case 1)
#
#
# min  f(X) = (1 - x)^2 + 100*(y - x^2)^2
#  X
#
# X = [x, y]^T
#
# s.t: --> (x - 1)^3 - y + 1 <= 0
#      --> x + y - 2 <= 0
#
# Bounds (or search domain):
# -1.5  <= x <= 1.5
# -0.5  <= y <= 2.5
#
# Global minimun:
# X*      = [1.0, 1.0]^T
# f(X*)   = 0
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
include("../../../01_-_Optimization_models/constrained_Rosenbrock_function_1.jl");

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by f and f_c respectivily

num_var, _ = size(bounds); # number of variables

vbounds = [-2.0*ones(num_var)  2.0*ones(num_var)]; # velocity limits

opt_arg = nothing;  # optional arguments
num_particles = 50; # number of particles
k_max = 200;        # maximum iteration
ω   = [0.9, 0.7];   # inertia weight
c1  = [1.0, 0.8];   # acceleration constants
c2  = [0.3, 0.5];

values_f, X_g_value, population_data, Swarm_data, fun_evals, const_evals = PSO_flyback(f, f_c, bounds,
                                                                                       opt_arg,
                                                                                       vbounds,
                                                                                       num_particles,
                                                                                       k_max, ω, c1, c2);

# solution
solution    = Swarm_data.X_global;
f_solution  = Swarm_data.f_X_global;

println("\n\n=================================================================")
println("Constrained Rosenbrock 1 function --> PSO flyback")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(0:k_max, values_f, xlabel="Iterations k",
             ylabel="f(x)", label="")) # plot
