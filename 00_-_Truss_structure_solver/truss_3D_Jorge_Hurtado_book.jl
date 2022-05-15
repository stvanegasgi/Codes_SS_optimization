# Structural model unknowns for 3D truss, example 2.2, Análisis matricial de
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
#   example 2.2, page 62
#
# -----------------------------------------------------------------------------

# node coordinates
#                     x        y          z
xcor =   Float64[     0.0     0.0       0.8216;
                     -2.5     0.0       0.6216;
                     -1.25   -2.165     0.6216;
                      1.25   -2.165     0.6216;
                      2.5     0.0       0.6216;
                      1.25    2.165     0.6216;
                     -1.25    2.165     0.6216;
                     -4.33   -2.5       0.0   ;
                      0.0    -5.0       0.0   ;
                      4.33   -2.5       0.0   ;
                      4.33    2.5       0.0   ;
                      0.0     5.0       0.0   ;
                     -4.33    2.5       0.0    ];

# elements information
#                N_start   N_end   Type_material
elements = Int64[   1        2           1;
                    1        3           1;
                    1        4           1;
                    1        5           1;
                    1        6           1;
                    1        7           1;
                    2        3           1;
                    3        4           1;
                    4        5           1;
                    5        6           1;
                    6        7           1;
                    2        7           1;
                    2        8           1;
                    3        8           1;
                    3        9           1;
                    4        9           1;
                    4       10           1;
                    5       10           1;
                    5       11           1;
                    6       11           1;
                    6       12           1;
                    7       12           1;
                    7       13           1;
                    2       13           1];

# materials properties
#                              Type   Areas        E
material_properties = Float64[  1    1.000e-04 2.058e+08];

# forces information (Direction X=1, Y=2, Z=3)
#                Node     Direction   Magnitude
forces = Float64[  1        3            -6;
                   2        3            -3;
                   3        3            -3;
                   4        3            -3;
                   5        3            -3;
                   6        3            -3;
                   7        3            -3];

# restricted degree of freedom (Direction X=1, Y=2, Z=3)
#                      Node   Direction
restricted_dof = Int64[  8       1;
                         9       1;
                        10       1;
                        11       1;
                        12       1;
                        13       1;
                         8       2;
                         9       2;
                        10       2;
                        11       2;
                        12       2;
                        13       2;
                         8       3;
                         9       3;
                        10       3;
                        11       3;
                        12       3;
                        13       3];
