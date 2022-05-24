# Structural model unknowns for 3D truss, load case 1 and 2 for truss with 72
# bars (case 1) with the PSOPC flyback method
#
#
# min  W(A) = Σ(ρ_i*L_i*A_i), for i=1, 2, ..., 72
#  A
#
# A = [A1, A2, ..., A72]^T
#        
# s.t:
# Stress limits: (compresive -) and (tensile +)
# Tensile limits:     all bars 25 [ksi]
# Compressive limits: all bar -25 [ksi]
#
# Limits displacement:
#-----------------------------------------
#           Displacement limits [in]
#   Node         x          y         z
#-----------------------------------------
#    17       +-0.25     +-0.25       -
#    18       +-0.25     +-0.25       -
#    19       +-0.25     +-0.25       -
#    20       +-0.25     +-0.25       -
#
# Bounds:
# --> (case 1) 0.1  [in^2] <= A_i <= 5.0 [in^2] for i=1, 2, ..., 72
#
# Number of restrictions: 320
# * Displacement: 2 (+-) * 2 Direction * 4 Nodes * 2 Cases  = 32
# * Stress: 2 (Compressive and Tensile) * 72 Bars * 2 Cases = 288
#
# keep in mind the units of unknows
# Areas:     [in^2]
# Forces:    [kips = 1_000 lb_f]
# Distances: [in]
# Stress:    [ksi = 1_000psi = 1_000 lb_f/in^2]
# Weight:    [lbf]
# Density:   [lb/in^3]
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
# *(1)  Li LJ, Huang ZB, Liu F, Wu QH. A heuristic particle swarm optimizer for
#       optimization of pin connected structures.
#       Comput Struct 2007;85(7–8):340–9
#
# -----------------------------------------------------------------------------

# ============================ packages =======================================

using Plots # for plots
ENV["GKSwstype"] = "100";

# ======================= PSOPC functions and model ===========================

include("../../../01_-_Optimization_models/truss_72_bars_3D_bound_1.jl"); # model
include("./PSOPC_flyback_functions.jl");                                  # load the PSOPC flyback functions

# ===================== data of optimization problem ==========================

# the objective and constraint function are defines by W and f_c respectivily

# variable bounds
bounds = truss_model_72_bars_bounds_1.variables_bounds;

num_var, _ = size(bounds); # number of variables

vmax = abs(5.0-0.1);
vbounds = [-vmax*ones(num_var)  vmax*ones(num_var)]; # velocity limits

opt_arg = truss_model_72_bars_bounds_1;  # optional arguments
num_particles = 50;                      # number of particles
k_max = 3000;                            # maximum iteration
ω   =  [0.9, 0.4];                       # inertia weight
c1  =  [0.8, 0.8];                       # acceleration constants
c2  =  [0.8, 0.8];
c3  =  [0.2, 0.4];

values_f, X_g_value, population_data, Swarm_data, fun_evals, const_evals = PSOPC_flyback(W, f_c, bounds,
                                                                                         opt_arg,
                                                                                         vbounds,
                                                                                         num_particles,
                                                                                         k_max, ω, c1, c2, c3);

# solution
solution    = Swarm_data.X_global;
f_solution  = Swarm_data.f_X_global;

println("\n\n=================================================================")
println("Truss 72 bars bounds 1 --> PSOPC flyback")
println("X* = $solution\n")
println("f(X*) = $f_solution\n")
println("Constraints (gi(X*) <= 0) = $(f_c(solution, opt_arg)) \n")
println("Objective function evaluations: $(fun_evals) \n")
println("Constraint function evaluations: $(const_evals) \n")

display(plot(0:k_max, values_f, xlabel="Iterations k",
             ylabel="f(x)", label="")) # plot
