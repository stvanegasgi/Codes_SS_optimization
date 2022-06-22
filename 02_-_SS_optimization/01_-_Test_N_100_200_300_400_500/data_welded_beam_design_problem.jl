# Welded beam design problem with the Subset Simulation (SS) optimization data
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

using JLD2;          # save variables
using Statistics;    # basic statistics functionality
using XLSX;          # save data in Excel

# ============================ data simulations ===============================

N_test = [100, 200, 300, 400, 500]; # samples to test

# variables to export to CSV
output_funtion          = zeros(4, length(N_test));
output_levels           = zeros(4, length(N_test));
output_evals_function   = zeros(4, length(N_test));
output_evals_constraint = zeros(4, length(N_test));
variables_mean          = zeros(4, length(N_test));
variables_std           = zeros(4, length(N_test));

for N in N_test

    file = "welded_beam_design_problem_N_"*string(N)*".jld2"
    @load file X_OPTIMAL_WBDP F_OPTIMAL_WBDP EVALS_F_OPTIMAL_WBDP EVALS_CONST_OPTIMAL_WBDP NUM_K_LEVELS_WBDP

    println("\nWelded beam design problem N $(N)")
    println("==========================================")
    println("Objective function: ")
    println("Best: ", minimum(F_OPTIMAL_WBDP))
    println("Mean: ", mean(F_OPTIMAL_WBDP))
    println("Worst: ", maximum(F_OPTIMAL_WBDP))
    println("Std: ", std(F_OPTIMAL_WBDP))

    output_funtion[:, Int64(N/100)] = [minimum(F_OPTIMAL_WBDP), mean(F_OPTIMAL_WBDP), maximum(F_OPTIMAL_WBDP), std(F_OPTIMAL_WBDP)];

    println("\nNumber of levels: ")
    println("Best: ", minimum(NUM_K_LEVELS_WBDP))
    println("Mean: ", mean(NUM_K_LEVELS_WBDP))
    println("Worst: ", maximum(NUM_K_LEVELS_WBDP))
    println("Std: ", std(NUM_K_LEVELS_WBDP))

    output_levels[:, Int64(N/100)] = [minimum(NUM_K_LEVELS_WBDP), mean(NUM_K_LEVELS_WBDP), maximum(NUM_K_LEVELS_WBDP), std(NUM_K_LEVELS_WBDP)];

    println("\nNumber of objective function evals: ")
    println("Best: ", minimum(EVALS_F_OPTIMAL_WBDP))
    println("Mean: ", mean(EVALS_F_OPTIMAL_WBDP))
    println("Worst: ", maximum(EVALS_F_OPTIMAL_WBDP))
    println("Std: ", std(EVALS_F_OPTIMAL_WBDP))

    output_evals_function[:, Int64(N/100)] = [minimum(EVALS_F_OPTIMAL_WBDP), mean(EVALS_F_OPTIMAL_WBDP), maximum(EVALS_F_OPTIMAL_WBDP), std(EVALS_F_OPTIMAL_WBDP)];

    println("\nNumber of constraint function evals: ")
    println("Best: ", minimum(EVALS_CONST_OPTIMAL_WBDP))
    println("Mean: ", mean(EVALS_CONST_OPTIMAL_WBDP))
    println("Worst: ", maximum(EVALS_CONST_OPTIMAL_WBDP))
    println("Std: ", std(EVALS_CONST_OPTIMAL_WBDP))

    output_evals_constraint[:, Int64(N/100)] = [minimum(EVALS_CONST_OPTIMAL_WBDP), mean(EVALS_CONST_OPTIMAL_WBDP), maximum(EVALS_CONST_OPTIMAL_WBDP), std(EVALS_CONST_OPTIMAL_WBDP)];

#   variables results
    variables_mean[:, Int64(N/100)]  = mean(X_OPTIMAL_WBDP, dims=2);
    variables_std[:, Int64(N/100)]   = std(X_OPTIMAL_WBDP, dims=2);
end

XLSX.openxlsx("welded_beam_design_problem.xlsx", mode="w") do xf

#   created the names sheets
    XLSX.rename!(xf[1], "Funtion_values");
    for sheetname in ["Level_values", "Evals_function", "Evals_constraint", "Variables_mean", "Variables_std"]
      XLSX.addsheet!(xf, sheetname);
    end

#   fill with the data
    for sheet in 1:4
      xf[sheet]["B1"] = "N = 100"; xf[sheet]["C1"] = "N = 200";
      xf[sheet]["D1"] = "N = 300"; xf[sheet]["E1"] = "N = 400";
      xf[sheet]["F1"] = "N = 500"; xf[sheet]["A1"] = "Index";
      xf[sheet]["A2"] = "Best";    xf[sheet]["A3"] = "Mean";
      xf[sheet]["A4"] = "Worst";   xf[sheet]["A5"] = "Std";
    end

    for sheet in 5:6
      xf[sheet]["B1"] = "N = 100"; xf[sheet]["C1"] = "N = 200";
      xf[sheet]["D1"] = "N = 300"; xf[sheet]["E1"] = "N = 400";
      xf[sheet]["F1"] = "N = 500"; xf[sheet]["A1"] = "Variables";
      xf[sheet]["A2"] = "h";       xf[sheet]["A3"] = "l";
      xf[sheet]["A4"] = "t";       xf[sheet]["A5"] = "b";
    end

    xf[1]["B2:F5"] = output_funtion;
    xf[2]["B2:F5"] = output_levels;
    xf[3]["B2:F5"] = output_evals_function;
    xf[4]["B2:F5"] = output_evals_constraint;
    xf[5]["B2:F5"] = variables_mean;
    xf[6]["B2:F5"] = variables_std;
end
