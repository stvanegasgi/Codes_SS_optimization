# Unconstrained Matyas function in two variables
#
#
# min  f(X) = 0.26*(x^2 + y^2) - 0.48xy
#  X
#
# X = [x, y]^T
#
# Bounds (or search domain):
# -10.0 <= x <= 10.0
# -10.0 <= y <= 10.0
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

# objective function (Matyas function)
function f(X::Array{Float64, 1}, opti::Any=nothing)

    x = X[1];   y = X[2];
    return 0.26*(x^2 + y^2) - 0.48*x*y
end;

# constraints function
function f_c(X::Array{Float64, 1}, opti::Any=nothing)

    return [-100] # unconstrained (always satisfy g(X) < 0)
end;

# variable bounds
#          min    max
bounds = [-10.0   10.0;  # --> x
          -10.0   10.0]; # --> y
