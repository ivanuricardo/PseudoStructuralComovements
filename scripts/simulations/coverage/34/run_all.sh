#!/bin/bash

julia -t 14 gamma_cov100.jl
julia -t 14 delta_cov100.jl
julia -t 14 gamma_cov250.jl
julia -t 14 delta_cov250.jl
