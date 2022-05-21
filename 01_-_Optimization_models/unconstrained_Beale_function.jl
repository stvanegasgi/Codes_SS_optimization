# Unconstrained Beale function in two variables
#
#
# min  f(X) = (1.5 - x + xy)^2 + (2.25 - x + xy^2)^2 + (2.625 - x + xy^3)^2
#  X
#
# X = [x, y]^T
#
# Bounds (or search domain):
# -4.5 <= x <= 4.5
# -4.5 <= y <= 4.5
#
# Global minimun:
# X*      = [3, 0.5]^T
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

# objective function (Beale function)
function f(X::Array{Float64, 1}, opti::Any=nothing)

    x = X[1];   y = X[2];
    return (1.5 - x + x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625 - x + x*y^3)^2
end;

# constraints function
function f_c(X::Array{Float64, 1}, opti::Any=nothing)

    return [-100] # unconstrained (always satisfy g(X) < 0)
end;

# variable bounds
#          min   max
bounds = [-4.5   4.5;  # --> x
          -4.5   4.5]; # --> y

# X* (optimal)
#             x    y
X_optimal = [3.0  0.5];
