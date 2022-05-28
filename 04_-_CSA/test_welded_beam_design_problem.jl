# Welded beam design problem with the crow search algorithm (CSA).
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

include("../01_-_Optimization_models/welded_beam_design_problem.jl"); # model
# load the CSA functions
include("./CSA_functions.jl");

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by f and f_c respectivily

# variable bounds
bounds = [0.1       2.0     # ----> x1
          0.1      10.0     # ----> x2
          0.1      10.0     # ----> x3
          0.1       2.0];   # ----> x4

X_min = bounds[:, 1]; # lower bounds
X_max = bounds[:, 2]; # upper bounds

opt_arg = nothing; # optional arguments
num_crows = 50;    # number of crows
k_max = 800;       # maximum iteration
fl = 2.0;          # flight length
AP = 0.1;          # awareness probability

x_optimal, f_x_optimal, fun_evals, const_evals = CSA(f, f_c, X_min, X_max, num_crows, k_max, fl, AP, opt_arg);

# solution
solution   = x_optimal[:, end];
f_solution = f_x_optimal[end];

println("\n\n=================================================================")
println("Welded beam design problem --> CSA")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(0:k_max, f_x_optimal, xlabel="Iterations k",
             ylabel="f(x)", label="")) # plot
