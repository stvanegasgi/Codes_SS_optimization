# Solved unconstrained Matyas function in two variables with the crow search algorithm (CSA)
#
#
# min  f(X) = 0.26*(x^2 + y^2) - 0.48xy
#  X
#
# X = [x, y]^T
#
# Bounds (or search domain):
# -10.0 <= x <= 10.0
# -10.0 <= y <= 10.0
#
# Global minimun:
# X*      = [0.0, 0.0]^T
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
include("../01_-_Optimization_models/unconstrained_Matyas_function.jl");

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by f and f_c respectivily

# variable bounds
#          min    max
bounds = [-10.0   10.0;  # --> x
          -10.0   10.0]; # --> y

X_min = bounds[:, 1]; # lower bounds
X_max = bounds[:, 2]; # upper bounds

opt_arg = nothing; # optional arguments
num_crows = 50;    # number of crows
k_max = 100;       # maximum iteration
fl = 2.0;          # flight length
AP = 0.1;          # awareness probability

x_optimal, f_x_optimal, fun_evals, const_evals = CSA(f, f_c, X_min, X_max, num_crows, k_max, fl, AP, opt_arg);

# solution
solution   = x_optimal[:, end];
f_solution = f_x_optimal[end];

println("\n\n=================================================================")
println("Unconstrained Matyas function --> CSA")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(0:k_max, f_x_optimal, xlabel="Iterations k",
            ylabel="f(x)", label="")) # plot