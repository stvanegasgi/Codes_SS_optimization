#!/bin/bash

echo "Runing all the case..."
julia test_welded_beam_design_problem.jl
julia test_3D_truss_25_bars.jl
julia test_3D_truss_72_bars_bound_1.jl
julia test_3D_truss_72_bars_bound_2.jl
julia test_2D_truss_200_bars.jl
echo "ok..."
