# Unconstrained Bukin function in two variables
#
#
# min  f(X) = 100*sqrt(abs(y - 0.01x^2)) + 0.01*abs(x + 10)
#  X
#
# X = [x, y]^T
#
# Bounds (or search domain):
# -15.0 <= x <= -5.0
#  -3.0 <= y <=  3.0
#
# Global minimun:
# X*      = [-10.0, 1.0]^T
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

# objective function (Bukin function)
function f(X::Array{Float64, 1}, opti::Any=nothing)

    x = X[1];   y = X[2];
    return 100*sqrt(abs(y - 0.01*x^2)) + 0.01*abs(x + 10)
end;

# constraints function
function f_c(X::Array{Float64, 1}, opti::Any=nothing)

    return [-100] # unconstrained (always satisfy g(X) < 0)
end;

# variable bounds
#          min     max
bounds = [-15.0   -5.0;  # --> x
          -3.0     3.0]; # --> y
