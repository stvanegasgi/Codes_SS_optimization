# Welded beam design problem with the Subset Simulation (SS) optimization
#
#
# min  f(X) = 1.10471*x1^2*x2 + 0.04811*x3*x4*(14.0 + x2)
#  X
#
# X = [x1, x2, x3, x4]^T = [h, l, t, b]^T
#
# s.t:
# ---> g1(X) = τ(X) - τ_max <= 0
# ---> g2(X) = σ(X) - σ_max <= 0
# ---> g3(X) = x1 - x4 <= 0
# ---> g4(X) = 0.10471*x1^2 + 0.04811*x3*x4*(14.0 + x2) - 5.0 <= 0
# ---> g5(X) = 0.125 - x1 <= 0
# ---> g6(X) = δ(X) - δ_max <= 0
# ---> g7(X) = P - Pc(X) <= 0
#
# where:
# --> τ(X) = sqrt(τ_p^2 + (2*τ_p*τ_pp*x2)/(2*R) + τ_pp^2);
# --> τ_p = P/(sqrt(2)*x1*x2);
# --> τ_pp = (M*R)/J;
# --> M = P*(L + x2/2);
# --> R = sqrt(x2^2/4 + ((x1 + x3)/2)^2);
# --> J = 2*(sqrt(2)*x1*x2*(x2^2/12 + ((x1 + x3)/2)^2));
# --> σ(X) = (6*P*L)/(x4*x3^2);
# --> δ(X) = (4*P*L^3)/(E*x3^2*x4);
# --> Pc(X) = ((4.013*E*sqrt((x3^2*x4^6)/36))/L^2)*(1-(x3/(2*L))*sqrt(E/(4*G)));
# --> P = 6_000 [lbf]; L = 14 [in]; E = 30*10^6 [psi]; G = 12*10^6 [psi];
# --> τ_max = 13_600 [psi]; σ_max = 30_000 [psi]; δ_max = 0.25 [in];
#
# Bounds:
# 0.1 <= x1 <=  2.0
# 0.1 <= x2 <= 10.0
# 0.1 <= x3 <= 10.0
# 0.1 <= x4 <=  2.0
#
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

# ======================= SS functions ========================================

include("../../01_-_Optimization_models/welded_beam_design_problem.jl"); # model
# load the SS functions
include("../ss_optimization.jl");

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by f and f_c respectivily

# variable bounds
bounds = [0.1       2.0     # ----> x1
          0.1      10.0     # ----> x2
          0.1      10.0     # ----> x3
          0.1       2.0];   # ----> x4

ε = 1e-6;          # convergence criterion
k_max = 3000;      # stop criterion
opt_arg = nothing; # optional argument
dimension = length(bounds[:, 1]); # dimension of the problem

executions = 30;   # executions

# preapred variables for saved (WBDP = welded beam design problem)
X_OPTIMAL_WBDP           = zeros(dimension, executions);
F_OPTIMAL_WBDP           = zeros(1, executions);
EVALS_F_OPTIMAL_WBDP     = zeros(1, executions);
EVALS_CONST_OPTIMAL_WBDP = zeros(1, executions);
NUM_K_LEVELS_WBDP        = zeros(1, executions);

N_test = [100, 200, 300, 400, 500]; # samples to test

for N in N_test

    println("\nWelded beam design problem N = $(N)")
    println("===============================================")

#   runs the problem
    for test_i = 1:executions

        x_optimal, f_x_optimal, samples_k_level, f_samples_k_level, hk_k_level, Fconk_k_level, fun_evals, const_evals = ss_optimization(f,
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
        X_OPTIMAL_WBDP[:, test_i]        = deepcopy(solution);
        F_OPTIMAL_WBDP[test_i]           = deepcopy(f_solution);
        EVALS_F_OPTIMAL_WBDP[test_i]     = deepcopy(fun_evals);
        EVALS_CONST_OPTIMAL_WBDP[test_i] = deepcopy(const_evals);
        NUM_K_LEVELS_WBDP[test_i]        = deepcopy(length(f_x_optimal));
    end

#   save the variables
    file = "welded_beam_design_problem_N_"*string(N)*".jld2"
    @save file X_OPTIMAL_WBDP F_OPTIMAL_WBDP EVALS_F_OPTIMAL_WBDP EVALS_CONST_OPTIMAL_WBDP NUM_K_LEVELS_WBDP
    println("The simulation welded beam design problem N = $(N) is OK...")
end
