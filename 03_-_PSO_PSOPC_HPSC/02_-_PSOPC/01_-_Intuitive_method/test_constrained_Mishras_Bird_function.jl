# Mishra's Bird function with the PSOPC flyback method
#
#
# min  f(X) = sin(y)e^((1 - cos(x))^2) + cos(x)e^((1 - sin(y))^2) + (x - y)^2
#  X
#
# X = [x, y]^T
#
# s.t: (x + 5)^2 + (y + 5)^2 - 25 < 0
#
# Bounds:
# -10  <= x <= 0
# -6.5 <= y <= 0
#
# Global minimun:
# X*      = [-3.1302468, -1.5821422]^T
# f(X*)   = -106.7645367
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

# ======================= PSOPC functions =====================================

# load the PSOPC flyback functions
include("./PSOPC_flyback_functions.jl")

# load the functions
include("../../../01_-_Optimization_models/constrained_Mishras_Bird_function.jl");

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by f and f_c respectivily

# variable bounds
bounds = [-10.0      0.0     # ----> x
           -6.5      0.0];   # ----> y

num_var, _ = size(bounds); # number of variables

vbounds = [-10.0*ones(num_var)   10.0*ones(num_var)]; # velocity limits

opt_arg = nothing;  # optional arguments
num_particles = 30; # number of particles
k_max = 200;        # maximum iteration
ω   = [0.7, 0.5];   # inertia weight
c1  = [2.0, 2.0];   # acceleration constants
c2  = [2.0, 2.0];
c3  = [0.3, 0.3];

values_f, X_g_value, population_data, Swarm_data = PSOPC_flyback(f, f_c, bounds,
                                                                 opt_arg,
                                                                 vbounds,
                                                                 num_particles,
                                                                 k_max, ω, c1, c2, c3);

# solution
solution    = Swarm_data.X_global;
f_solution  = Swarm_data.f_X_global;

println("\n\n=================================================================")
println("Constrained Mishras Bird function --> PSOPC flyback")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(0:k_max, values_f, xlabel="Iterations k",
             ylabel="f(x)", label="")) # plot
