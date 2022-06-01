# Mishra's Bird function with with the Subset Simulation (SS) optimization
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
# DATE:    June 2022
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

# ============================ packages =======================================

using Plots # for plots
ENV["GKSwstype"] = "100";

# ======================= SS functions ========================================

# load the SS functions
include("./ss_optimization.jl")

# load the functions
include("../01_-_Optimization_models/constrained_Mishras_Bird_function.jl");

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by f and f_c respectivily

# variable bounds
bounds = [-10.0      0.0     # ----> x
           -6.5      0.0];   # ----> y

N = 100;           # number of samples
ε = 1e-5;          # convergence criterion
opt_arg = nothing; # optional argument

x_optimal, f_x_optimal, samples_k_level, f_samples_k_level, hk_k_level, Fconk_k_level, fun_evals, const_evals = ss_optimization(f,
                                                                                                                                f_c,
                                                                                                                                N,
                                                                                                                                bounds,
                                                                                                                                opt_arg,
                                                                                                                                ε);

# solution
solution    = x_optimal[:, end];
f_solution  = f_x_optimal[end];

println("\n\n=================================================================")
println("Constrained Mishras Bird function --> SS")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(1:length(f_x_optimal), f_x_optimal, xlabel="Iterations k",
            ylabel="f(x)", label="")) # plot
