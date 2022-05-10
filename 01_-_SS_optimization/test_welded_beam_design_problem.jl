# Welded beam design problem with the SS optimization method. Reference [1]
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
# The reference [1] has some errors in the constraints and objective function,
# see the reference [2]
#
# =============================================================================
# ------- April 2021
# -----------------------------------------------------------------------------
# by 
# ------- Steven Vanegas Giraldo -------> stvanegasgi@unal.edu.co
# -----------------------------------------------------------------------------
# Universidad Nacional de Colombia - Sede Manizales
# =============================================================================
#
# References:
# *(1) Li, H.-S., Au, S.-K (2010). Desing optimization using subset simulation
#      algorithm. Structural Safety, 32(6), 384-392.
#
# *(2) Coello Coello CA. Use of a self-adaptive penalty approach for engineering
#      optimization problems. Comput Ind 2000;41(2):13–127.

# -----------------------------------------------------------------------------


# ======================= SS optimization functions ===========================

include("./ss_optimization.jl"); # load the SS optimization functions

# ===================== data of optimization problem ==========================

# objective function (Welded beam cost function)
function f(X, opti)

    x1 = X[1];  # h
    x2 = X[2];  # l
    x3 = X[3];  # t
    x4 = X[4];  # b
    return -(1.10471*x1^2*x2 + 0.04811*x3*x4*(14.0 + x2))
end

# constraints function
function f_c(X, opti)

    x1 = X[1];  # h
    x2 = X[2];  # l
    x3 = X[3];  # t
    x4 = X[4];  # b

    P = 6_000;        # [lbf]
    L = 14;           # [in]
    E = 30e6;         # [psi]
    G = 12e6;         # [psi]
    τ_max = 13_600;   # [psi]
    σ_max = 30_000;   # [psi]
    δ_max = 0.25;     # [in]

    M    = P*(L + x2/2);
    R    = sqrt(x2^2/4 + ((x1 + x3)/2)^2);
    J    = 2*(sqrt(2)*x1*x2*(x2^2/12 + ((x1 + x3)/2)^2));
    τ_p  = P/(sqrt(2)*x1*x2);
    τ_pp = (M*R)/J;
    τ    = sqrt(τ_p^2 + (2*τ_p*τ_pp*x2)/(2*R) + τ_pp^2);
    σ    = (6*P*L)/(x4*x3^2);
    δ    = (4*P*L^3)/(E*x3^3*x4);
    Pc   = ((4.013*E*sqrt((x3^2*x4^6)/36))/L^2)*(1-(x3/(2*L))*sqrt(E/(4*G)));

    g1 = τ - τ_max;                                         # g1(X) <= 0
    g2 = σ - σ_max;                                         # g2(X) <= 0
    g3 = x1 - x4;                                           # g3(X) <= 0
    g4 = 0.10471*x1^2 + 0.04811*x3*x4*(14.0 + x2) - 5.0;    # g4(X) <= 0
    g5 = 0.125 - x1;                                        # g5(X) <= 0
    g6 = δ - δ_max;                                         # g6(X) <= 0
    g7 = P - Pc;                                            # g7(X) <= 0

    return [g1, g2, g3, g4, g5, g6, g7]
end

# variable bounds
bounds = [0.1       2.0     # ----> x1
          0.1      10.0     # ----> x2
          0.1      10.0     # ----> x3
          0.1       2.0];   # ----> x4

opt_arg = nothing;
ε = 1e-4; # convergence criterion (10^-5  -  10^-7)

x_opti, h_opti = ss_optimization(f, f_c, 100, bounds, opt_arg, ε);

println(x_opti[:, end])
println(h_opti)

