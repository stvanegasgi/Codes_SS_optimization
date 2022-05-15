# Unconstrained Goldstein-Price function in two variables
#
#
# min  f(X) = (1 + (x + y + 1)^2 (19 - 14x + 3x^2 - 14y + 6xy + 3y^2) (30 + (2x - 3y)^2 (18 - 32x + 12x^2 + 48y - 36xy + 27y^2)
#  X
#
# X = [x, y]^T
#
# Bounds (or search domain):
# -2.0 <= x <= 2.0
# -2.0 <= y <= 2.0
#
# Global minimun:
# X*      = [0.0, -1.0]^T
# f(X*)   = 3
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

# objective function (Goldstein-Price function)
function f(X::Array{Float64, 1}, opti::Any=nothing)

    x = X[1];   y = X[2];
    term_1 = (1 + (x + y + 1)^2 * (19 - 14*x +  3*x^2 - 14*y +  6*x*y +  3*y^2));
    term_2 = (30 + (2*x - 3*y)^2 * (18 - 32*x + 12*x^2 + 48*y - 36*x*y + 27*y^2));
    return term_1 * term_2
end;

# constraints function
function f_c(X::Array{Float64, 1}, opti::Any=nothing)

    return [-100] # unconstrained (always satisfy g(X) < 0)
end;

# variable bounds
#          min   max
bounds = [-2.0   2.0;  # --> x
          -2.0   2.0]; # --> y
