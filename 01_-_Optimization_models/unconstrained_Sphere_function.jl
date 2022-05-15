# Unconstrained Sphere function in multiple variables
#
#             n
# min  f(X) = Σ xi^2
#  X         i=1
#
# X = [x1, ..., xn]^T
#
# Bounds (or search domain):
# -∞ <= xi <= ∞
#
# Global minimun:
# X*      = [0.0, ..., 0.0]^T
# f(X*)   = 0
#
# =============================================================================
# DATE:    November 2021
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

# objective function (Sphere function)
function f(X::Array{Float64, 1}, opti::Any=nothing)

    return sum(X.^2)
end;

# constraints function
function f_c(X::Array{Float64, 1}, opti::Any=nothing)

    return [-100] # unconstrained (always satisfy g(X) < 0)
end;

# variable bounds
#           min               max
bounds = [-1.0e20*ones(n) 1.0e20*ones(n)];   # -∞ <= xi <= ∞
