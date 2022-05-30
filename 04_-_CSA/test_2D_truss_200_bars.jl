# Structural model unknowns for 3D truss, load case 1, 2 and 3 for truss with
# 200 bars with the crow search algorithm (CSA)
#
#
# min  W(A) = Σ(ρ_i*L_i*A_i), for i=1, 2, ..., 200
#  A
#
# A = [A1, A2, ..., A200]^T
#
# s.t:
# Stress limits: (compresive -) and (tensile +)
# Tensile limits:     all bars 10 [ksi]
# Compressive limits: all bar -10 [ksi]
#
# Limits displacement:
# No restricted
#
# Bounds:
# --> 0.1  [in^2] <= A_i <= 20.0 [in^2] for i=1, 2, ..., 200
#
# Number of restrictions: 1200
# * Displacement: = 0
# * Stress: 2 (Compressive and Tensile) * 200 Bars * 3 load case = 1200
#
# keep in mind the units of unknows
# Areas:     [in^2]
# Forces:    [kips = 1_000 lb_f]
# Distances: [in]
# Stress:    [ksi = 1_000psi = 1_000 lb_f/in^2]
# Weight:    [lbf]
# Density:   [lb/in^3]
#
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

# ======================= CSA functions and model =============================

include("../01_-_Optimization_models/truss_200_bars_2D.jl"); # model
include("./CSA_functions.jl");                               # load the CSA functions

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by W and f_c respectivily

# variable bounds
bounds = truss_model_200_bars.variables_bounds;

X_min = bounds[:, 1]; # lower bounds
X_max = bounds[:, 2]; # upper bounds

opt_arg = truss_model_200_bars; # optional arguments
num_crows = 60;                 # number of crows
k_max = 2500;                   # maximum iteration
fl = 2.5;                       # flight length
AP = 0.1;                       # awareness probability

x_optimal, f_x_optimal, fun_evals, const_evals = CSA(W, f_c, X_min, X_max, num_crows, k_max, fl, AP, opt_arg);

# solution
solution   = x_optimal[:, end];
f_solution = f_x_optimal[end];

println("\n\n=================================================================")
println("Truss 200 bars --> CSA")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(0:k_max, f_x_optimal, xlabel="Iterations k",
             ylabel="f(x)", label="")) # plot
