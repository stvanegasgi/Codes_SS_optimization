# Structural model unknowns for 2D truss, example 11.3, Análisis estructural
# Jairo Uribe Escamilla, segunda edición
#
#
# keep in mind the units of unknows
# Areas:     [cm^2]
# Forces:    [tn_f]
# Distances: [cm]
# Stress:    [tn_f/cm^2]
#
# =============================================================================
# ------- November 2021
# -----------------------------------------------------------------------------
# by 
# ------- Steven Vanegas Giraldo -------> stvanegasgi@unal.edu.co
# -----------------------------------------------------------------------------
# Universidad Nacional de Colombia - Sede Manizales
# =============================================================================
#
# References:
# * Análisis estructural - Jairo Uribe Escamilla - 2 edición
#   example 11.3, page 437
#
# -----------------------------------------------------------------------------

# node coordinates
#                     x       y
xcor =   Float64[    0.0     0.0;
                   800.0     0.0;
                   400.0   300.0;
                   400.0     0.0];

# elements information
#                N_start   N_end   Type_material
elements = Int64[   1        3         1;
                    1        4         2;
                    3        2         3;
                    4        2         2;
                    3        4         4];

# materials properties
#                              Type   Areas     E
material_properties = Float64[  1     100     2040;
                                2      40     2040;
                                3     150     2040;
                                4      30     2040];

# forces information (Direction X=1, Y=2, Z=3)
#                Node     Direction   Magnitude
forces = Float64[ 4           2         -20;
                  3           1           4;
                  3           2           3];

# restricted degree of freedom (Direction X=1, Y=2, Z=3)
#                      Node   Direction
restricted_dof = Int64[  1       1;
                         1       2;
                         2       2];
