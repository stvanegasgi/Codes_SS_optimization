# Unconstrained Rosembrock function in multiple variables
#
#            n-1
# min  f(X) = Σ (100*(x(i+1) - xi^2)^2 + (1 - xi)^2)
#  X         i=1
#
# X = [x1, ..., xn]^T
#
# Bounds (or search domain):
# -∞ <= xi <= ∞
#
# Global minimun:
# X*      = [1.0, ..., 1.0]^T
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

# objective function (Rosembrock function)
function f(X::Array{Float64, 1}, opti::Any=nothing)

    n = length(X);
    fun = 0;
    for i in 1:n-1
        fun += 100*(X[i+1] - X[i]^2)^2 + (1 - X[i])^2;
    end
    return fun
end;

# constraints function
function f_c(X::Array{Float64, 1}, opti::Any=nothing)

    return [-100] # unconstrained (always satisfy g(X) < 0)
end;

# variable bounds
#              min              max
bounds = [-1.0e20*ones(n)  1.0e20*ones(n)];   # -∞ <= xi <= ∞
