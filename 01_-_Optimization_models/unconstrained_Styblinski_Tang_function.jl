# Unconstrained Styblinski-Tang function in two variables
#
#                     n
# min  f(X) = (1/2) * Σ (xi^4 - 16xi^2 + 5xi)
#  X                 i=1
#
# X = [x, y]^T
#
# Bounds (or search domain):
# -5.0 <= xi <= 5.0
#
# Global minimun:
# X*      = [-2.903534, ..., -2.903534]^T
# -39.16617n <= f(X*) <= -39.16616n
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

# ============================== varibales ====================================

n = 2; # number of variables (change for multiple variables)

# ===================== data of optimization problem ==========================

# objective function (Styblinski-Tang function)
function f(X::Array{Float64, 1}, opti::Any=nothing)

    return (1/2) * sum(X.^4 - 16*X.^2 + 5*X)
end;

# constraints function
function f_c(X::Array{Float64, 1}, opti::Any=nothing)

    return [-100] # unconstrained (always satisfy g(X) < 0)
end;

# variable bounds
#               min        max
bounds = [-5.0*ones(n) 5.0*ones(n)];
