using LinearAlgebra
using Format
using Base.Threads
using Statistics
using BenchmarkTools
using Revise
includet("HillClimbZeroOne.jl")
using .HillClimbZeroOne


function __main__()

	####################
	# Variables & Main #
	####################

	my_config = ClimbConfig(
		n_size = 22,			# n x n matrix
		c_climbs = 30000,		# max climbs per trial
		t_trials = 100,		# total trials (random starts) || Distributed if parallel
		tabu_max = 30,			# max memory for tabu list
		parallel = false,
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

    HillClimbZeroOne.simulate_climb(my_config)

end

Base.invokelatest(__main__)