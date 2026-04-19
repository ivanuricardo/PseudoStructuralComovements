#!/bin/bash

# julia -t 14 comparison_delta100.jl
# julia -t 14 comparison_delta250.jl
print "running gamma"
julia -t 14 comparison_gamma250.jl
julia -t 14 comparison_gamma100.jl
