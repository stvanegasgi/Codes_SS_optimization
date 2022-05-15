# Unconstrained Rastrigin function in multiple variables
#
#                   n
# min  f(X) = A*n + Σ (xi^2 - A*cos(2πxi)), with A = 10
#  X               i=1
#
# X = [x1, ..., xn]^T
#
# Bounds (or search domain):
# -5.12 <= xi <= 5.12
#
# Global minimun:
# X*      = [0.0, ..., 0.0]^T
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

# ============================== varibales ====================================

n = 2; # number of variables (change for multiple variables)

# ===================== data of optimization problem ==========================

# objective function (Rastrigin function)
function f(X::Array{Float64, 1}, opti::Any=nothing)

    n = length(X);  A = 10;
    return A*n + sum(X.^2 - A*cos.(2*π*X))
end;

# constraints function
function f_c(X::Array{Float64, 1}, opti::Any=nothing)

    return [-100] # unconstrained (always satisfy g(X) < 0)
end;

# variable bounds
#               min         max
bounds = [-5.12*ones(n) 5.12*ones(n)];
