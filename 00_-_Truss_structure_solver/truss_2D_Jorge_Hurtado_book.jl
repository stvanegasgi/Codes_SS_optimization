# Structural model unknowns for 2D truss, example 2.1, Análisis matricial de
# estructuras Jorge Eduardo Hurtado (Curso con Matlab)
#
#
# keep in mind the units of unknows
# Areas:     [m^2]
# Forces:    [kN]
# Distances: [m]
# Stress:    [kN/m^2]
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
# * Análisis matricial de estructuras - Jorge Eduardo Hurtado
#   example 2.1, page 49
#
# -----------------------------------------------------------------------------

# node coordinates
#                     x       y
xcor =   Float64[    0.0     0.0;
                     4.0     0.0;
                     4.0     3.0;
                     8.0     0.0;
                     8.0     3.0;
                    12.0     0.0];

# elements information
#                N_start   N_end   Type_material
elements = Int64[   1        2          1;
                    1        3          1;
                    2        3          1;
                    2        4          1;
                    3        4          1;
                    3        5          1;
                    4        5          1;
                    4        6          1;
                    5        6          1];

# materials properties
#                              Type   Areas     E
material_properties = Float64[   1    5e-03   2e+08];

# forces information (Direction X=1, Y=2, Z=3)
#                Node     Direction   Magnitude
forces = Float64[  3           1          4;
                   3           2        -20;
                   5           2        -20];

# restricted degree of freedom (Direction X=1, Y=2, Z=3)
#                      Node   Direction
restricted_dof = Int64[  1        1;
                         1        2;
                         6        1;
                         6        2];
