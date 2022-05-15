# Structural model unknowns for 3D truss, example 11.4, Análisis estructural
# Jairo Uribe Escamilla, segunda edición
#
#
# keep in mind the units of unknows
# Areas:     [mm^2]
# Forces:    [kN]
# Distances: [mm]
# Stress:    [kN/mm^2]
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
# * Análisis estructural - Jairo Uribe Escamilla - 2 edición
#   example 11.4, page 444
#
# -----------------------------------------------------------------------------

# node coordinates
#                     x         y         z
xcor =   Float64[   2250.0    6000.0    4800.0;
                    3750.0    6000.0    2400.0;
                    5250.0    6000.0    4800.0;
                       0.0       0.0    6000.0;
                    3750.0       0.0       0.0;
                    7500.0       0.0    6000.0];

# elements information
#                N_start   N_end   Type_material
elements = Int64[   1        2          1;
                    1        3          1;
                    1        4          2;
                    1        6          3;
                    2        3          1;
                    2        4          3;
                    2        5          2;
                    3        5          3;
                    3        6          2;
                    4        5          4;
                    4        6          4;
                    5        6          4];

# materials properties
#                              Type   Areas     E
material_properties = Float64[   1     2000    210;
                                 2     4000    210;
                                 3     5000    210;
                                 4     1000    210];

# forces information (Direction X=1, Y=2, Z=3)
#                Node     Direction   Magnitude
forces = Float64[  1          1         100;
                   1          2         150;
                   1          3        -120;
                   2          1          50;
                   2          2         -30;
                   2          3        -100;
                   3          1         -40;
                   3          2         -20;
                   3          3         -60];

# restricted degree of freedom (Direction X=1, Y=2, Z=3)
#                      Node   Direction
restricted_dof = Int64[  4        1;
                         4        2;
                         4        3;
                         5        1;
                         5        2;
                         5        3;
                         6        1;
                         6        2;
                         6        3];
