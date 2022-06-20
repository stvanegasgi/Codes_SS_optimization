# Structural model unknowns for 3D truss, load case 1 and 2 for truss with 72
# bars (case 1) with the Subset Simulation (SS) optimization
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

# ======================= SS functions and model ==============================

include("../01_-_Optimization_models/truss_72_bars_3D_bound_1.jl"); # model
include("./ss_optimization.jl");                                    # load the SS

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by W and f_c respectivily

# variable bounds
bounds = truss_model_72_bars_bounds_1.variables_bounds;

N = 100;           # number of samples
ε = 1e-3;          # convergence criterion
k_max = 1000;      # stop criterion
opt_arg = truss_model_72_bars_bounds_1; # optional arguments

x_optimal, f_x_optimal, samples_k_level, f_samples_k_level, hk_k_level, Fconk_k_level, fun_evals, const_evals = ss_optimization(W,
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
println("Truss 72 bars bounds 1 --> SS")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(1:length(f_x_optimal), f_x_optimal, xlabel="Iterations k",
             ylabel="f(x)", label="")) # plot
