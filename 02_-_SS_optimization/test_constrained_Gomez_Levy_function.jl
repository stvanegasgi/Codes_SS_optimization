# Solve Gomez and Levy function (modified) with the Subset Simulation (SS) optimization
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
include("./ss_optimization.jl");

# load the functions
include("../01_-_Optimization_models/constrained_Gomez_Levy_function.jl");

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by f and f_c respectivily

# variable bounds
bounds = [-1.0   0.75;    # ----> x
          -1.0   1.0 ];   # ----> y

N = 100;           # number of samples
ε = 1e-5;          # convergence criterion
k_max = 3000;      # stop criterion
opt_arg = nothing; # optional argument

x_optimal, f_x_optimal, samples_k_level, f_samples_k_level, hk_k_level, Fconk_k_level, fun_evals, const_evals = ss_optimization(f,
                                                                                                                                f_c,
                                                                                                                                N,
                                                                                                                                bounds,
                                                                                                                                opt_arg,
                                                                                                                                ε,
                                                                                                                                k_max);

# solution
solution    = x_optimal[:, end];
f_solution  = f_x_optimal[end];

println("\n\n=================================================================")
println("Constrained Gomez Levy function --> SS")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(1:length(f_x_optimal), f_x_optimal, xlabel="Iterations k",
             ylabel="f(x)", label="")) # plot
