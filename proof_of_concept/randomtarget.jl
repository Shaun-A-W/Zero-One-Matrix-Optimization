using LinearAlgebra
include("sw-01matrix-search.jl")


max_n = 20
max_determinants = Dict(
    1 => 1,
    2 => 1,
    3 => 2,
    4 => 3,
    5 => 5,
    6 => 9,
    7 => 32,
    8 => 56,
    9 => 144,
    10 => 320,
    11 => 1458,
    12 => 3645,
    13 => 9477,
    14 => 25515,
    15 => 131072,
    16 => 327680,
    17 => 1114112,
    18 => 3411968,
    19 => 19531250,
    20 => 56640625,
    21 => 195312500,
    22 => 967396000,
)

# n_size = 15				# n x n matrix
c_climbs = 1_000_0		# max climbs per trials
t_trials = 1_000_0     	# total trials (random starts)
tabu_max = 50			# max memory for tabu list
tabu_list = Tuple{Int, Int}[]	# empty tabu list, do not touch

move_counts = Int64[]
for n_size in 1:20
    move_counts = Int64[]
    fail_count = 0
    println("Simulating random targetting on $n_size by $n_size matrix.")
    for i in 1:100
        target = rand(1:max_determinants[n_size])
        # print("Target: $target  //  ")
        best_target, best_prox, end_move = main_loop_target(n_size, c_climbs, t_trials, tabu_max, tabu_list, target)
        if isnothing(end_move)
            fail_count += 1
            continue
        else
            push!(move_counts, end_move)
        end
    end
    avg_c = sum(move_counts) / length(move_counts)
    println("Average moves: $avg_c")
    println("Fails: $fail_count")
end

