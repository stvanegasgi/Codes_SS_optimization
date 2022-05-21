# Unconstrained Ackley function in two variables
#
#
# min  f(X) = -20*e^(-0.2*sqrt(0.5(x^2 + y^2))) - e^(0.5(cos(2πx) + cos(2πy))) + e + 20
#  X
#
# X = [x, y]^T
#
# Bounds (or search domain):
# -5 <= x <= 5
# -5 <= y <= 5
#
# Global minimun:
# X*      = [0.0, 0.0]^T
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
# *(1) https://en.wikipedia.org/wiki/Test_functions_for_optimization
#
# -----------------------------------------------------------------------------

# ===================== data of optimization problem ==========================

# objective function (Ackley function)
function f(X::Array{Float64, 1}, opti::Any=nothing)

    x = X[1];   y = X[2];
    return -20*exp(-0.2*sqrt(0.5*(x^2 + y^2))) - exp(0.5*(cos(2*π*x) + cos(2*π*y))) + exp(1) + 20
end;

# constraints function
function f_c(X::Array{Float64, 1}, opti::Any=nothing)

    return [-100] # unconstrained (always satisfy g(X) < 0)
end;

# variable bounds
#          min      max
bounds = [-5.0      5.0;    # ----> x
          -5.0      5.0];   # ----> y

# X* (optimal)
#             x    y
X_optimal = [0.0  0.0];
