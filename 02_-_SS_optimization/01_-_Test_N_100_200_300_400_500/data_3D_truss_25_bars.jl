# Structural model unknowns for 3D truss, load case 1 and 2 for truss with 25
# bars with the Subset Simulation (SS) optimization data
#
#
# min  W(A) = Σ(ρ_i*L_i*A_i), for i=1, 2, ..., 25
#  A
#
# A = [A1, A2, ..., A25]^T
#        
# s.t:
# Stress limits:
#-------------------------------------------------------------------------
#    Variables         Compressive stress [ksi]      Tensile stress [ksi]
#-------------------------------------------------------------------------
#   1       A1                 35.092                      40.0
#   2    A2-A5                 11.590                      40.0
#   3    A6-A9                 17.305                      40.0
#   4  A10-A11                 35.092                      40.0
#   5  A12-A13                 35.092                      40.0
#   6  A14-A17                  6.759                      40.0
#   7  A18-A21                  6.959                      40.0
#   8  A22-A25                 11.082                      40.0
#
# All limits displacement in three direction: +- 0.35 [in]
#
# Bounds:
# 0.01 [in^2] <= Ai <= 3.5 [in], for i=1, 2, ..., 25
#
# Number of restrictions: 172
# * Displacement: 2 (+-) * 3 Direction * 6 Nodes * 2 Cases  = 72
# * Stress: 2 (Compressive and Tensile) * 25 Bars * 2 Cases = 100
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
variables_mean          = zeros(8, length(N_test));
variables_std           = zeros(8, length(N_test));

for N in N_test

    file = "truss_25_bars_N_"*string(N)*".jld2"
    @load file X_OPTIMAL_TRUSS_25 F_OPTIMAL_TRUSS_25 EVALS_F_OPTIMAL_TRUSS_25 EVALS_CONST_OPTIMAL_TRUSS_25 NUM_K_LEVELS_TRUSS_25

    println("\nTruss 25 bars problem N $(N)")
    println("==========================================")
    println("Objective function: ")
    println("Best: ", minimum(F_OPTIMAL_TRUSS_25))
    println("Mean: ", mean(F_OPTIMAL_TRUSS_25))
    println("Worst: ", maximum(F_OPTIMAL_TRUSS_25))
    println("Std: ", std(F_OPTIMAL_TRUSS_25))

    output_funtion[:, Int64(N/100)] = [minimum(F_OPTIMAL_TRUSS_25), mean(F_OPTIMAL_TRUSS_25), maximum(F_OPTIMAL_TRUSS_25), std(F_OPTIMAL_TRUSS_25)];

    println("\nNumber of levels: ")
    println("Best: ", minimum(NUM_K_LEVELS_TRUSS_25))
    println("Mean: ", mean(NUM_K_LEVELS_TRUSS_25))
    println("Worst: ", maximum(NUM_K_LEVELS_TRUSS_25))
    println("Std: ", std(NUM_K_LEVELS_TRUSS_25))

    output_levels[:, Int64(N/100)] = [minimum(NUM_K_LEVELS_TRUSS_25), mean(NUM_K_LEVELS_TRUSS_25), maximum(NUM_K_LEVELS_TRUSS_25), std(NUM_K_LEVELS_TRUSS_25)];

    println("\nNumber of objective function evals: ")
    println("Best: ", minimum(EVALS_F_OPTIMAL_TRUSS_25))
    println("Mean: ", mean(EVALS_F_OPTIMAL_TRUSS_25))
    println("Worst: ", maximum(EVALS_F_OPTIMAL_TRUSS_25))
    println("Std: ", std(EVALS_F_OPTIMAL_TRUSS_25))

    output_evals_function[:, Int64(N/100)] = [minimum(EVALS_F_OPTIMAL_TRUSS_25), mean(EVALS_F_OPTIMAL_TRUSS_25), maximum(EVALS_F_OPTIMAL_TRUSS_25), std(EVALS_F_OPTIMAL_TRUSS_25)];

    println("\nNumber of constraint function evals: ")
    println("Best: ", minimum(EVALS_CONST_OPTIMAL_TRUSS_25))
    println("Mean: ", mean(EVALS_CONST_OPTIMAL_TRUSS_25))
    println("Worst: ", maximum(EVALS_CONST_OPTIMAL_TRUSS_25))
    println("Std: ", std(EVALS_CONST_OPTIMAL_TRUSS_25))

    output_evals_constraint[:, Int64(N/100)] = [minimum(EVALS_CONST_OPTIMAL_TRUSS_25), mean(EVALS_CONST_OPTIMAL_TRUSS_25), maximum(EVALS_CONST_OPTIMAL_TRUSS_25), std(EVALS_CONST_OPTIMAL_TRUSS_25)];

#   variables results
    variables_mean[:, Int64(N/100)]  = mean(X_OPTIMAL_TRUSS_25, dims=2);
    variables_std[:, Int64(N/100)]   = std(X_OPTIMAL_TRUSS_25, dims=2);
end

XLSX.openxlsx("truss_25_bars_problem.xlsx", mode="w") do xf

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
      xf[sheet]["A2"] = "A(1)";    xf[sheet]["A3"] = "A(2)";
      xf[sheet]["A4"] = "A(3)";    xf[sheet]["A5"] = "A(4)";
      xf[sheet]["A6"] = "A(5)";    xf[sheet]["A7"] = "A(6)";
      xf[sheet]["A8"] = "A(7)";    xf[sheet]["A9"] = "A(8)";
    end

    xf[1]["B2:F5"] = output_funtion;
    xf[2]["B2:F5"] = output_levels;
    xf[3]["B2:F5"] = output_evals_function;
    xf[4]["B2:F5"] = output_evals_constraint;
    xf[5]["B2:F9"] = variables_mean;
    xf[6]["B2:F9"] = variables_std;
end
