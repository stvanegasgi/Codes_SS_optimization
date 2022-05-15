# Structural model unknowns for 3D truss, load case 1 and 2 for truss with 72
# bars (case 1 and case 2) with subset simulation optimization, bounds case 1
#
#
# min  W(A) = Σ(ρ_i*L_i*A_i), for i=1, 2, ..., 72
#  A
#
# A = [A1, A2, ..., A72]^T
#        
# s.t:
# Stress limits: (compresive -) and (tensile +)
# Tensile limits: all bars 25 [ksi]
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
#
# =============================================================================
# ------- April 2022
# -----------------------------------------------------------------------------
# by 
# ------- Steven Vanegas Giraldo -------> stvanegasgi@unal.edu.co
# -----------------------------------------------------------------------------
# Universidad Nacional de Colombia - Sede Manizales
# =============================================================================
#
# References:
# *(1) Li, H.-S., Au, S.-K (2010). Desing optimization using subset simulation
#      algorithm. Structural Safety, 32(6), 384-392.
#
# *(2) Haftka RT, Gurdal Z. Elements of structural optimization. 3rd
#      ed. Dordrecht: Kluwer Academic Publishers; 1992. (pag-245)
#
# -----------------------------------------------------------------------------


# ======================= SS functions and model ==============================

include("./load_case_1_2_3D_truss_72_bars_bounds_1.jl"); # load the variables of model  
include("./structure_solver_functions.jl");              # function structural solver
include("./ss_optimization.jl");                         # load the SS optimization functions

# ===================== data of optimization problem ==========================

# objective function (weight of structure)
function W(A, model)

    L_i           = model.length_bars; # the length of each bar
    Type_material = model.info_elements[:, 3]; # type of material
    ρ_i           = model.info_material_propieties[Type_material, 4]; # density 
    A_i           = A[Type_material]; # areas for each element

    return -sum(ρ_i .* L_i .* A_i); # for min problem
end

# constraints function
function f_c(A, model)

#   material propieties with A to test
    material_properties_load_1_2 = deepcopy(model.info_material_propieties[:, 1:3]);
    material_properties_load_1_2[:, 2] = A;

#   solution FEM load case 1
    δ1, _, _, σ1, _, _ = solver_truss(model.coordinates,
                                      model.info_elements,
                                      material_properties_load_1_2,
                                      model.load_cases[1],
                                      model.restricted_dof);

#   stress constraints for the load case 1
    σ1_T = σ1 - model.axial_sigma_restricted[:, 2];   # σ - σmax <= 0
    σ1_C = model.axial_sigma_restricted[:, 3]  - σ1;  # σmin - σ <= 0

#   displacement constraint
    nodes           = model.constraint_displacement[:, 1];
    direction       = model.constraint_displacement[:, 2];
    constraint_dof  = Int64.(3*nodes .- 3 .+ direction);

#   positive and negative displacement
#   δ - δmax <= 0
    δ1_p = δ1[constraint_dof] - model.constraint_displacement[:, 3];
#   δmin - δ <= 0
    δ1_n = model.constraint_displacement[:, 4] - δ1[constraint_dof];

#   solution FEM load case 2
    δ2, _, _, σ2, _, _ = solver_truss(model.coordinates,
                                      model.info_elements,
                                      material_properties_load_1_2,
                                      model.load_cases[2],
                                      model.restricted_dof);

#   stress constraints for the load case 2
    σ2_T = σ2 - model.axial_sigma_restricted[:, 2];   # σ - σmax <= 0
    σ2_C = model.axial_sigma_restricted[:, 3]  - σ2;  # σmin - σ <= 0

#   positive and negative displacement
#   δ - δmax <= 0
    δ2_p = δ2[constraint_dof] - model.constraint_displacement[:, 3];
#   δmin - δ <= 0
    δ2_n = model.constraint_displacement[:, 4] - δ2[constraint_dof];

    return [σ1_T; σ1_C; δ1_p; δ1_n; σ2_T; σ2_C; δ2_p; δ2_n]
end

# variable bounds
bounds = truss_model_72_bars_bounds_1.variables_bounds;
opt_arg = truss_model_72_bars_bounds_1;  # optional arguments
ε = 1e-6; # convergence criterion (10^-5  -  10^-7)

x_opti, h_opti = ss_optimization(W, f_c, 300, bounds, opt_arg, ε);

println(x_opti[:, end])
println(h_opti)
