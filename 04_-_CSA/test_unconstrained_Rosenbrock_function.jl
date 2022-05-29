# Solved unconstrained Rosembrock function in multiple variables with the crow search algorithm (CSA)
# method
#
#            n-1
# min  f(X) = Σ (100*(x(i+1) - xi^2)^2 + (1 - xi)^2)
#  X         i=1
#
# X = [x1, ..., xn]^T
#
# Bounds (or search domain):
# -∞ <= xi <= ∞
#
# Global minimun:
# X*      = [1.0, ..., 1.0]^T
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
include("../01_-_Optimization_models/unconstrained_Rosenbrock_function.jl");

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by f and f_c respectivily

# variable bounds
#          min      max
bounds = [-2.0      2.0;
          -2.0      3.0];   # -∞ <= xi <= ∞ (for plotting and examples)

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
println("Unconstrained Rosembrock function --> CSA")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(0:k_max, f_x_optimal, xlabel="Iterations k",
            ylabel="f(x)", label="")) # plot
