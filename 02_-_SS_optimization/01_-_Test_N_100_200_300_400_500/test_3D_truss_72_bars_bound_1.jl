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

using JLD2; # save variables

# ======================= SS functions and model ==============================

include("../../01_-_Optimization_models/truss_72_bars_3D_bound_1.jl"); # model
include("../ss_optimization.jl");                                      # load the SS

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by W and f_c respectivily

# variable bounds
bounds = truss_model_72_bars_bounds_1.variables_bounds;

ε = 1e-4;          # convergence criterion
k_max = 2000;      # stop criterion
opt_arg = truss_model_72_bars_bounds_1; # optional arguments
dimension = length(bounds[:, 1]); # dimension of the problem

executions = 30;   # executions

# preapred variables for saved
X_OPTIMAL_TRUSS_72_1           = zeros(dimension, executions);
F_OPTIMAL_TRUSS_72_1           = zeros(1, executions);
EVALS_F_OPTIMAL_TRUSS_72_1     = zeros(1, executions);
EVALS_CONST_OPTIMAL_TRUSS_72_1 = zeros(1, executions);
NUM_K_LEVELS_TRUSS_72_1        = zeros(1, executions);

N_test = [100, 200, 300, 400, 500]; # samples to test

for N in N_test

    println("\nTruss 72 bars case 1 problem N = $(N)")
    println("===============================================")

#   runs the problem
    for test_i = 1:executions

        x_optimal, f_x_optimal, samples_k_level, f_samples_k_level, hk_k_level, Fconk_k_level, fun_evals, const_evals = ss_optimization(W,
                                                                                                                                        f_c,
                                                                                                                                        N,
                                                                                                                                        bounds,
                                                                                                                                        opt_arg,
                                                                                                                                        ε,
                                                                                                                                        k_max);

#       solution
        solution    = x_optimal[:, end];
        f_solution  = f_x_optimal[end];

#       save data for the simulation
        X_OPTIMAL_TRUSS_72_1[:, test_i]        = deepcopy(solution);
        F_OPTIMAL_TRUSS_72_1[test_i]           = deepcopy(f_solution);
        EVALS_F_OPTIMAL_TRUSS_72_1[test_i]     = deepcopy(fun_evals);
        EVALS_CONST_OPTIMAL_TRUSS_72_1[test_i] = deepcopy(const_evals);
        NUM_K_LEVELS_TRUSS_72_1[test_i]        = deepcopy(length(f_x_optimal));
    end

#   save the variables
    file = "truss_72_bars_case_1_N_"*string(N)*".jld2"
    @save file X_OPTIMAL_TRUSS_72_1 F_OPTIMAL_TRUSS_72_1 EVALS_F_OPTIMAL_TRUSS_72_1 EVALS_CONST_OPTIMAL_TRUSS_72_1 NUM_K_LEVELS_TRUSS_72_1
    println("The simulation truss 72 bars case 1 problem N = $(N) is OK...")
end
