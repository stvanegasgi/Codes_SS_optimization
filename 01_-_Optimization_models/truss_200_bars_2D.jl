# Structural model unknowns for 3D truss, load case 1, 2 and 3 for truss with
# 200 bars
#
#             n
# min  W(A) = Σ(ρ_i*L_i*A_i), for i=1, 2, ..., 200
#  A         i=1
#
# A = [A1, A2, ..., A200]^T
#
# s.t:
# Stress limits: (compresive -) and (tensile +)
# Tensile limits:     all bars 10 [ksi]
# Compressive limits: all bar -10 [ksi]
#
# Limits displacement:
# No restricted
#
# Bounds:
# --> 0.1  [in^2] <= A_i <= 20.0 [in^2] for i=1, 2, ..., 200
#
# Number of restrictions: 1200
# * Displacement: = 0
# * Stress: 2 (Compressive and Tensile) * 200 Bars * 3 load case = 1200
#
# keep in mind the units of unknows
# Areas:     [in^2]
# Forces:    [kips = 1_000 lb_f]
# Distances: [in]
# Stress:    [ksi = 1_000psi = 1_000 lb_f/in^2]
# Weight:    [lbf]
# Density:   [lb/in^3]
#
#
# =============================================================================
# DATE:    December 2021
# WHO:     Steven Vanegas Giraldo
# EMAIL:   stvanegasgi@unal.edu.co
# -----------------------------------------------------------------------------
# Universidad Nacional de Colombia - Sede Manizales
# =============================================================================
#
# References:
# *(1) Awad, R. (2021, October). Sizing optimization of truss structures using
#      the political optimizer (PO) algorithm. In Structures
#      (Vol. 33, pp. 4871-4894). Elsevier.
#
# *(2) Degertekin, S. O., & Hayalioglu, M. S. (2013). Sizing truss structures 
#      using teaching-learning-based optimization.
#      Computers & Structures, 119, 177-188.
#
# -----------------------------------------------------------------------------

# ============== load the model data and truss solver =========================

include("./load_case_1_2_3_2D_truss_200_bars.jl"); # load the variables of model
include("../00_-_Truss_structure_solver/structure_solver_functions.jl"); # load the structural solver

# ===================== data of optimization problem ==========================

# objective function (weight of structure)
function W(A, model)

    L_i           = model.length_bars; # the length of each bar
    Type_material = model.info_elements[:, 3]; # type of material
    ρ_i           = model.info_material_propieties[Type_material, 4]; # density 
    A_i           = A[Type_material]; # areas for each element

    return sum(ρ_i .* L_i .* A_i);
end;

# constraints function
function f_c(A, model)

#   material propieties with A to test
    material_properties_load_1_2_3 = copy(model.info_material_propieties[:, 1:3]);
    material_properties_load_1_2_3[:, 2] = A;

#   solution FEM load case 1
    _, _, _, σ1, _, _ = solver_truss(model.coordinates,
                                     model.info_elements,
                                     material_properties_load_1_2_3,
                                     model.load_cases[1],
                                     model.restricted_dof);

#   stress constraints for the load case 1
    σ1_T = σ1 - model.axial_sigma_restricted[:, 2];   # σ - σmax <= 0
    σ1_C = model.axial_sigma_restricted[:, 3] - σ1;  # σmin - σ <= 0

#   solution FEM load case 2
    _, _, _, σ2, _, _ = solver_truss(model.coordinates,
                                      model.info_elements,
                                      material_properties_load_1_2_3,
                                      model.load_cases[2],
                                      model.restricted_dof);

#   stress constraints for the load case 2
    σ2_T = σ2 - model.axial_sigma_restricted[:, 2];   # σ - σmax <= 0
    σ2_C = model.axial_sigma_restricted[:, 3] - σ2;  # σmin - σ <= 0

#   solution FEM load case 3
    _, _, _, σ3, _, _ = solver_truss(model.coordinates,
                                     model.info_elements,
                                     material_properties_load_1_2_3,
                                     model.load_cases[3],
                                     model.restricted_dof);

#   stress constraints for the load case 1
    σ3_T = σ3 - model.axial_sigma_restricted[:, 2];   # σ - σmax <= 0
    σ3_C = model.axial_sigma_restricted[:, 3] - σ3;  # σmin - σ <= 0

    return [σ1_T; σ1_C; σ2_T; σ2_C; σ3_T; σ3_C]
end;

# variable bounds
bounds = truss_model_200_bars.variables_bounds;
