# Solve Townsend function (modified) in two variables with the Subset Simulation (SS) optimization
#
#
#
# min  f(X) = -(cos((x - 0.1) * y))^2 - x*sin(3x + y)
#  X
#
# X = [x, y]^T
#
# s.t: x^2 + y^2 -(2cos(t) - 0.5cos(2t) - 0.25cos(3t) - 0.125cos(4t))^2 - (2sin(t))^2 < 0
# t = atan2(x,y)
#
# Bounds:
# -2.25 <= x <= 2.25
# -2.5  <= y <= 1.75
#
# Global minimun:
# X*     = [2.0052938, 1.1944509]^T
# f(X*)  = -2.0239884
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
include("../01_-_Optimization_models/constrained_Townsend_function.jl");

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by f and f_c respectivily

# variable bounds
bounds = [-2.25   2.25;    # ----> x
          -2.5    1.75];   # ----> y

N = 200;           # number of samples
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
println("Constrained Townsend function --> SS")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(1:length(f_x_optimal), f_x_optimal, xlabel="Iterations k",
            ylabel="f(x)", label="")) # plot
