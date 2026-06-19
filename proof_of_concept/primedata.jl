using LinearAlgebra
using Primes
using Base.Threads
include("sw-01matrix-search.jl")


# p = 2693
# position = length(primes(p))
# print(position)

n_size = 12				# n x n matrix
c_climbs = 1_000_0		# max climbs per trials
t_trials = 1_000_0     	# total trials (random starts)
tabu_max = 20			# max memory for tabu list
tabu_list = Tuple{Int, Int}[]	# empty tabu list, do not touch

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


for n in 1:20
    println("N: $n, Total Primes in (Theoretical) Bounds: $(length(primes(max_determinants[n])))")
end

first_prime_fail = Dict{Int, Tuple{String, String}}()

# Lock prevents threads from printing over each other in the console
print_lock = ReentrantLock()
# Lock protects dictionary from data races
dict_lock = ReentrantLock()

start_time = time()
println("Available Threads: $(Sys.CPU_THREADS)")
println("Threads in use: $(Threads.nthreads())")
# Assign each n x n case to available cores (depends on IDE, VSCode in settings.json, etc.)
@threads for n in 1:max_n

    lock(print_lock) do
        println("Core $(threadid()) is starting prime analysis for $n x $n case...")
    end
    
    max_det = max_determinants[n]
    n_primes = primes(max_det)
    
    # Run the inner loop on assigned core
    for (pos, target) in enumerate(n_primes)
        best_target, best_prox, end_move = main_loop_target(n, c_climbs, t_trials, tabu_max, tabu_list, target)
        
        if best_prox != 0
            lock(print_lock) do
                println("Core $(threadid()) found first failure for $n at Prime: $target (Pos: $pos).")
            end
            
            lock(dict_lock) do
                first_prime_fail[n] = ("Prime: $target", "Position: $pos")
            end
            
            break # Exit after first fail
        end
    end
end

#     for target in n_primes
#         best_target, best_prox, end_move = main_loop_target(n, c_climbs, t_trials, tabu_max, tabu_list, target)
#         if best_prox != 0
#             println("Failure at prime determinant: $target.")
#             pos = length(primes(target))
#             println("Prime position = $pos.")
#             first_prime_fail[n] = ("Prime: $target", "Position: $pos")
#             break
#         end
#     end
# end

println(first_prime_fail)
println("Test concluded.")
println("Time taken: $(time() - start_time) seconds.")