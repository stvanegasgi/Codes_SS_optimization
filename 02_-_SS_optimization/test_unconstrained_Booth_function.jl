# Solved unconstrained Booth function in two variables with the Subset Simulation (SS) optimization
#
#
# min  f(X) = (x + 2y - 7)^2 + (2x + y - 5)^2
#  X
#
# X = [x, y]^T
#
# Bounds (or search domain):
# -10 <= x <= 10
# -10 <= y <= 10
#
# Global minimun:
# X*      = [1.0, 3.0]^T
# f(X*)   = 0
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
include("../01_-_Optimization_models/unconstrained_Booth_function.jl");

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by f and f_c respectivily

# variable bounds
#          min         max
bounds = [-10.0      10.0;    # ----> x
          -10.0      10.0];   # ----> y

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
println("Unconstrained Booth function --> SS")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(1:length(f_x_optimal), f_x_optimal, xlabel="Iterations k",
            ylabel="f(x)", label="")) # plot
