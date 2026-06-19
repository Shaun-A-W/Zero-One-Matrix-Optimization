using LinearAlgebra
using Primes
include("sw-01matrix-search.jl")

p_list = primes(3645)
fails = Set([])

n_size = 12				# n x n matrix
c_climbs = 1_000		# max climbs per trials
t_trials = 1_000		# total trials (random starts)
tabu_max = 20			# max memory for tabu list
tabu_list = Tuple{Int, Int}[]	# empty tabu list, do not touch

println(length(p_list))
println("Starting primes target search on $n_size x $n_size case...")
start_time = time()
for target in p_list
    # start_time = time()
	# println("Starting new search for $target on a $n_size x $n_size Matrix.")
	best_target, best_prox = main_loop_target(n_size, c_climbs, t_trials, tabu_max, tabu_list, target)
    if best_prox != 0
        # println("Determinant $target failed to find.")
        push!(fails, target)
    # else
    #     println("Determinant $target successfully found.")
    end
end
end_time = time()
println("Time taken: ", cfmt("%.2f", end_time - start_time), " seconds")
println("Compiled fails: ", fails)


# println("Now search for nxn prime fails on n+1 x n+1 case...")
# println("(12x12 max = 3645)")
# for target in fails
#     best_target, best_prox = main_loop_target(n_size+1, c_climbs, t_trials, tabu_max, tabu_list, target)
#     if best_prox != 0
#         println("Determinant $target failed to find.")
#     end
# end
# println("Test end.")