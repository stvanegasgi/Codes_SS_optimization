# Structural model unknowns for 3D truss, load case 1 and 2 for truss with 72
# bars (case 2)
#
#             n
# min  W(A) = Σ(ρ_i*L_i*A_i), for i=1, 2, ..., 72
#  A         i=1
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
# --> (case 2) 0.01  [in^2] <= A_i <= 5.0 [in^2] for i=1, 2, ..., 72
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
# *(1) Li, H.-S., Au, S.-K (2010). Desing optimization using subset simulation
#      algorithm. Structural Safety, 32(6), 384-392.
#
# *(2) L. A. Schmit and H. Miura. Approximation concepts for efficient
#      structural synthesis, Tech. Report CR-2552, NASA, 1976.
#
# -----------------------------------------------------------------------------

# ============== load the model data and truss solver =========================

include("./load_case_1_2_3D_truss_72_bars_bounds_2.jl"); # load the variables of model
include("../00_-_Truss_structure_solver/structure_solver_functions.jl"); # load the structural solver

# ===================== data of optimization problem ==========================

# objective function (weight of structure)
function W(A, model)

#   lecture index
    TYPE_MAT_IND = 3;   DENSITY_IND = 4;

    L_i           = model.length_bars; # the length of each bar
    Type_material = model.info_elements[:, TYPE_MAT_IND]; # type of material
    ρ_i           = model.info_material_propieties[Type_material, DENSITY_IND]; # density 
    A_i           = A[Type_material]; # areas for each element

    return sum(ρ_i .* L_i .* A_i);
end;

# constraints function
function f_c(A, model)

#   lecture index
    AREA_IND = 2;   TENSILE_IND = 2;    COMPRESIVE_IND = 3;    DMAX_IND = 3;
    DMIN_IND = 4;

#   material propieties with A to test
    material_properties_load_1_2 = copy(model.info_material_propieties[:, 1:3]);
    material_properties_load_1_2[:, AREA_IND] = A;

#   solution FEM load case 1
    δ1, _, _, σ1, _, _ = solver_truss(model.coordinates,
                                      model.info_elements,
                                      material_properties_load_1_2,
                                      model.load_cases[1],
                                      model.restricted_dof);

#   stress constraints for the load case 1
    σ1_T = σ1 - model.axial_sigma_restricted[:, TENSILE_IND];     # σ - σmax <= 0
    σ1_C = model.axial_sigma_restricted[:, COMPRESIVE_IND]  - σ1; # σmin - σ <= 0

#   displacement constraint
    nodes           = model.constraint_displacement[:, 1];
    direction       = model.constraint_displacement[:, 2];
    constraint_dof  = Int64.(3*nodes .- 3 .+ direction);

#   positive and negative displacement
#   δ - δmax <= 0
    δ1_p = δ1[constraint_dof] - model.constraint_displacement[:, DMAX_IND];
#   δmin - δ <= 0
    δ1_n = model.constraint_displacement[:, DMIN_IND] - δ1[constraint_dof];

#   solution FEM load case 2
    δ2, _, _, σ2, _, _ = solver_truss(model.coordinates,
                                      model.info_elements,
                                      material_properties_load_1_2,
                                      model.load_cases[2],
                                      model.restricted_dof);

#   stress constraints for the load case 2
    σ2_T = σ2 - model.axial_sigma_restricted[:, TENSILE_IND];     # σ - σmax <= 0
    σ2_C = model.axial_sigma_restricted[:, COMPRESIVE_IND]  - σ2; # σmin - σ <= 0

#   positive and negative displacement
#   δ - δmax <= 0
    δ2_p = δ2[constraint_dof] - model.constraint_displacement[:, DMAX_IND];
#   δmin - δ <= 0
    δ2_n = model.constraint_displacement[:, DMIN_IND] - δ2[constraint_dof];

    return [σ1_T; σ1_C; δ1_p; δ1_n; σ2_T; σ2_C; δ2_p; δ2_n]
end;

# variable bounds
bounds = truss_model_72_bars_bounds_2.variables_bounds;
