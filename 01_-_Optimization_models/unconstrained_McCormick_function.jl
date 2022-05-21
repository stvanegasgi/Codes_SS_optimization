# Unconstrained McCormick function in two variables
#
#
# min  f(X) = sin(x + y) + (x - y)^2 - 1.5x + 2.5y + 1
#  X
#
# X = [x, y]^T
#
# Bounds (or search domain):
# -1.5 <= x <= 4.0
# -3.0 <= y <= 4.0
#
# Global minimun:
# X*      = [-0.54719, -1.54719]^T
# f(X*)   = -1.9133
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

# objective function (McCormick function)
function f(X::Array{Float64, 1}, opti::Any=nothing)

    x = X[1];   y = X[2];
    return sin(x + y) + (x - y)^2 - 1.5*x + 2.5*y + 1
end;

# constraints function
function f_c(X::Array{Float64, 1}, opti::Any=nothing)

    return [-100] # unconstrained (always satisfy g(X) < 0)
end;

# variable bounds
#          min    max
bounds = [-1.5   4.0;  # --> x
          -3.0   4.0]; # --> y

# X* (optimal)
#               x         y
X_optimal = [-0.54719  -1.54719];
