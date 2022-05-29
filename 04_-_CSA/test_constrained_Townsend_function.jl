# Solve Townsend function (modified) in two variables with the crow search algorithm (CSA)
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
# DATE:    May 2022
# WHO:     Steven Vanegas Giraldo
# EMAIL:   stvanegasgi@unal.edu.co
# -----------------------------------------------------------------------------
# Universidad Nacional de Colombia - Sede Manizales
# =============================================================================
#
# References:
# * Askarzadeh, A. (2016). A novel metaheuristic method for solving constrained
#   engineering optimization problems: crow search algorithm.
#   Computers & Structures, 169, 1-12.
#
# -----------------------------------------------------------------------------

# ============================ packages =======================================

using Plots # for plots
ENV["GKSwstype"] = "100";

# ======================= CSA functions =======================================

# load the CSA functions
include("./CSA_functions.jl");

# load the functions
include("../01_-_Optimization_models/constrained_Townsend_function.jl");

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by f and f_c respectivily

# variable bounds
bounds = [-2.25   2.25;    # ----> x
          -2.5    1.75];   # ----> y

X_min = bounds[:, 1]; # lower bounds
X_max = bounds[:, 2]; # upper bounds

opt_arg = nothing; # optional arguments
num_crows = 80;    # number of crows
k_max = 1000;      # maximum iteration
fl = 2.0;          # flight length
AP = 0.1;          # awareness probability

x_optimal, f_x_optimal, fun_evals, const_evals = CSA(f, f_c, X_min, X_max, num_crows, k_max, fl, AP, opt_arg);

# solution
solution   = x_optimal[:, end];
f_solution = f_x_optimal[end];

println("\n\n=================================================================")
println("Constrained Townsend function --> CSA")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(0:k_max, f_x_optimal, xlabel="Iterations k",
             ylabel="f(x)", label="")) # plot
