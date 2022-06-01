# Solved unconstrained McCormick function in two variables with the Subset Simulation (SS) optimization
#
#
# min  f(X) = sin(x + y) + (x - y)^2 - 1.5x + 2.5y + 1
#  X
#
# X = [x, y]^T
#
# Bounds (or search domain):
# -1.5 <= x <= 4.0
# -3.0 <= y <= 4.0
#
# Global minimun:
# X*      = [-0.54719, -1.54719]^T
# f(X*)   = -1.9133
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
include("../01_-_Optimization_models/unconstrained_McCormick_function.jl");

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by f and f_c respectivily

# variable bounds
#          min    max
bounds = [-1.5   4.0;  # --> x
          -3.0   4.0]; # --> y

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
println("Unconstrained McCormick function --> SS")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(1:length(f_x_optimal), f_x_optimal, xlabel="Iterations k",
            ylabel="f(x)", label="")) # plot
