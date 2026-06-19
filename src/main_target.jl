using LinearAlgebra
using Format
using Base.Threads
using Statistics
using BenchmarkTools
using Revise
includet("TargetClimbZeroOne.jl")
using .TargetClimbZeroOne


function __main__()

	####################
	# Variables & Main #
	####################

	my_config = TargetConfig(
		target = 2000,			# target determinant value
		n_size = 22,			# n x n matrix
		c_climbs = 30000,		# max climbs per trial
		t_trials = 10000,		# total trials (random starts) || Distributed if parallel
		tabu_max = 30,			# max memory for tabu list
		parallel = true,
		verbose = true,
		seed = nothing,
		thread_count = 0,
		output = true
	)

	# for seeding template 
	A =[1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;
		0  0  1  0  0  1  1  0  0  0  1  0  1  1  0  0  1  1  1  1  1  1;
		0  0  0  1  0  0  1  1  1  1  1  0  1  0  0  1  1  0  1  1  0  0;
		0  0  1  0  1  0  1  0  1  1  1  1  0  1  1  0  0  1  1  1  0  0;
		0  0  0  1  1  1  0  1  0  1  0  0  1  1  1  0  0  1  1  0  0  1;
		0  0  1  1  1  0  1  0  0  0  1  1  1  0  1  1  0  0  1  0  1  1;
		0  0  0  0  0  1  0  0  1  1  1  1  1  0  1  1  1  1  0  0  0  1;
		0  1  1  1  1  1  1  1  1  0  1  0  0  0  0  0  0  1  0  0  0  1;
		0  1  1  0  1  0  0  0  1  1  0  0  1  0  0  1  0  1  1  1  1  1;
		0  0  1  0  1  1  0  1  0  1  1  0  0  0  1  1  1  0  0  1  1  0;
		0  0  0  1  1  0  1  0  1  0  0  0  0  1  1  1  1  1  0  1  1  1;
		0  1  1  1  0  0  0  1  0  0  0  1  1  0  1  0  1  1  0  1  0  0;
		0  1  1  0  1  1  1  0  1  0  0  0  1  1  1  1  1  0  1  0  0  0;
		0  1  0  0  0  0  1  1  0  1  1  0  1  1  1  1  0  1  0  0  1  0;
		0  0  1  0  1  0  1  1  1  1  0  1  1  1  0  0  1  0  0  0  1  1;
		0  1  0  1  1  1  0  0  1  0  1  1  1  1  0  0  0  0  0  1  1  0;
		0  1  1  1  0  0  0  0  1  1  1  0  0  1  1  0  1  0  1  0  1  1;
		0  1  1  1  0  1  1  0  0  1  0  1  0  1  0  1  0  0  0  1  0  1;
		0  1  0  0  0  1  1  1  1  0  0  1  0  0  1  0  0  0  1  1  1  1;
		0  0  1  1  0  1  0  1  1  0  0  1  0  1  0  1  0  1  1  0  1  0;
		0  1  0  0  1  0  0  1  0  0  1  1  0  1  0  1  1  0  1  1  0  1;
		0  1  0  1  1  1  1  0  0  1  0  1  0  0  0  0  1  1  1  0  1  0]
	
	# println("d a ", cfmt("%.2f", det(A)))

	TargetClimbZeroOne.simulate_target(my_config)

end

Base.invokelatest(__main__)