# Structural model unknowns for 3D truss, load case 1 and 2 for truss with 72
# bars (case 1) with the crow search algorithm (CSA)
#
#
# min  W(A) = Σ(ρ_i*L_i*A_i), for i=1, 2, ..., 72
#  A
#
# A = [A1, A2, ..., A72]^T
#        
# s.t:
# Stress limits: (compresive -) and (tensile +)
# Tensile limits:     all bars 25 [ksi]
# Compressive limits: all bar -25 [ksi]
#
# Limits displacement:
#-----------------------------------------
#           Displacement limits [in]
#   Node         x          y         z
#-----------------------------------------
#    17       +-0.25     +-0.25       -
#    18       +-0.25     +-0.25       -
#    19       +-0.25     +-0.25       -
#    20       +-0.25     +-0.25       -
#
# Bounds:
# --> (case 1) 0.1  [in^2] <= A_i <= 5.0 [in^2] for i=1, 2, ..., 72
#
# Number of restrictions: 320
# * Displacement: 2 (+-) * 2 Direction * 4 Nodes * 2 Cases  = 32
# * Stress: 2 (Compressive and Tensile) * 72 Bars * 2 Cases = 288
#
# keep in mind the units of unknows
# Areas:     [in^2]
# Forces:    [kips = 1_000 lb_f]
# Distances: [in]
# Stress:    [ksi = 1_000psi = 1_000 lb_f/in^2]
# Weight:    [lbf]
# Density:   [lb/in^3]
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

include("../01_-_Optimization_models/truss_72_bars_3D_bound_1.jl"); # model
include("./CSA_functions.jl");                                      # load the CSA

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by W and f_c respectivily

# variable bounds
bounds = truss_model_72_bars_bounds_1.variables_bounds;

X_min = bounds[:, 1]; # lower bounds
X_max = bounds[:, 2]; # upper bounds

opt_arg = truss_model_72_bars_bounds_1; # optional arguments
num_crows = 60;                         # number of crows
k_max = 2500;                           # maximum iteration
fl = 2.5;                               # flight length
AP = 0.1;                               # awareness probability

x_optimal, f_x_optimal, fun_evals, const_evals = CSA(W, f_c, X_min, X_max, num_crows, k_max, fl, AP, opt_arg);

# solution
solution   = x_optimal[:, end];
f_solution = f_x_optimal[end];

println("\n\n=================================================================")
println("Truss 72 bars bounds 1 --> CSA")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(0:k_max, f_x_optimal, xlabel="Iterations k",
             ylabel="f(x)", label="")) # plot
