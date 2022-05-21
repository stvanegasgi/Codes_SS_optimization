# Unconstrained Booth function in two variables
#
#
# min  f(X) = (x + 2y - 7)^2 + (2x + y - 5)^2
#  X
#
# X = [x, y]^T
#
# Bounds (or search domain):
# -10 <= x <= 10
# -10 <= y <= 10
#
# Global minimun:
# X*      = [1.0, 3.0]^T
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

# objective function (Booth function)
function f(X::Array{Float64, 1}, opti::Any=nothing)

    x = X[1];   y = X[2];
    return (x + 2*y - 7)^2 + (2*x + y - 5)^2
end;

# constraints function
function f_c(X::Array{Float64, 1}, opti::Any=nothing)

    x = X[1];   y = X[2];
    return [-100] # unconstrained (always satisfy g(X) < 0)
end;

# variable bounds
#          min         max
bounds = [-10.0      10.0;    # ----> x
          -10.0      10.0];   # ----> y

# X* (optimal)
#             x     y
X_optimal = [1.0   3.0];
