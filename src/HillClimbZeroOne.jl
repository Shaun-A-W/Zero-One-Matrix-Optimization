module HillClimbZeroOne

# Shaun Worthen
# UTulsa Department of Mathematics
# Analysis of Zero-One Matrices
# Optimization / Search Algorithms

using LinearAlgebra
using Format
using Base.Threads
using Statistics
using Revise

export simulate_climb
export matrix_gradient
export ClimbConfig


Base.@kwdef struct ClimbConfig
    n_size::Int
    c_climbs::Int
    t_trials::Int
    tabu_max::Int
    seed::Union{Matrix{Int}, Nothing} = nothing
    verbose::Bool = false
    parallel::Bool = false
    thread_count::Int = 0
    output::Bool = true
end

struct Workspace{T<:Integer, F<:AbstractFloat}
    gradient::Matrix{T}
    A_precise::Matrix{F}
    A_inv::Matrix{F}
    I_mat::Matrix{F}
    tabu_mat::BitMatrix
    tabu_order_queue::Vector{Tuple{Int,Int}}
end

# A helper function to initialize the workspace cleanly
function initialize_workspace(n::Int, tabu_max::Int, ::Type{T}, ::Type{F}) where {T, F}
    return Workspace{T, F}(
        zeros(T, n, n),
        zeros(F, n, n),
        zeros(F, n, n),
        Matrix{F}(I, n, n),
        falses(n, n),
        Tuple{Int, Int}[]
    )
end


# Calculate gradient (change in determinant)
function matrix_gradient(A::Matrix{Int}, current_det::Real, ws::Workspace{T, F}, tabu_max::Int) where {T, F}
    (; gradient, A_precise, A_inv, I_mat, tabu_mat, tabu_order_queue) = ws
    n = size(A, 1)
	d = current_det

	# If singular, force invertibility
    # Inconsequential due to wanting high determinant
    while abs(d) < 1e-9
        shuffle_bits(A, n, 1, tabu_mat, tabu_order_queue, tabu_max)
        d = det(A)
    end

    # if abs(d) < 1e-9
        # grad = zeros(Int, n, n)
        # for i in 1:n
        #     for j in 1:n
        #         # Using @views avoids memory allocations during slicing
        #         sub_matrix = @views A[1:n .!= i, 1:n .!= j]
        #         sub_det = n < 23 ? Int64(round(det(sub_matrix))) : LinearAlgebra.det_bareiss(T.(sub_matrix))
        #         grad[i, j] = ((-1)^(i + j)) * Int64(sub_det)
        #     end
        # end
        # return grad
    
    # If invertible, use shortcut
    # Mind the memory optimization...

    A_precise .= A
    d_precise = T == BigInt ? BigFloat(d) : Float64(d)
    # in-place LU Factorization of A_precise
    lu_A = lu!(A_precise; check=false) 
    # reset A_inv to the Identity matrix
    A_inv .= I_mat 
    # solve A_precise * A_inv = I in-place
    ldiv!(lu_A, A_inv)
    # zero-allocation loop for the scalar multiplication and transpose
    for j in 1:n
        for i in 1:n
            gradient[i, j] = round(Int, d_precise * A_inv[j, i]) 
        end
    end

    return T(round(d))
end


function shuffle_bits(A::Matrix{Int}, n::Int, shakes::Int, tabu_mat::BitMatrix, tabu_order_queue::Vector{Tuple{Int,Int}}, tabu_max::Int)
    for _ in 1:shakes
        i, j = rand(1:n), rand(1:n)
        A[i, j] = 1 - A[i, j]

        # maintain set of moves in bool matrix
        tabu_mat[i, j] = true
        # maintain order of moves
        push!(tabu_order_queue, (i, j)) 

        # grab oldest move and discard
        if length(tabu_order_queue) > tabu_max
            old_i, old_j = popfirst!(tabu_order_queue)
            tabu_mat[old_i, old_j] = false # Free the oldest Tabu move
        end
    end
end


# Determine best move to increase determinant
function _best_change(A::Matrix{Int}, tabu_mat::BitMatrix, gradient::Matrix{T}, current_det::Real, best_det::Real) where {T<:Integer}
    n = size(A, 1)
    best_move = (row = 0, col = 0)
    best_delta = -Inf
    found_any = false

	# Note: Flipping 0 -> 1 adds grad
	#		Flipping 1 -> 0 subs grad
    
    for j in 1:n
        for i in 1:n
            # find change
            delta = A[i, j] == 0 ? gradient[i, j] : -gradient[i, j]
            
            # skip worse moves than current best
            if delta > best_delta
                # check Tabu memory
                is_tabu = tabu_mat[i, j] # Instantaneous lookup

                # check if new best or not Tabu
                if (current_det + delta) > best_det || !is_tabu
                    best_delta = delta
                    best_move = (row = i, col = j)
                    found_any = true
                end
            end
        end
    end
    
    return best_move, best_delta, found_any
end


# Primary internal loop used per trial
# Performed once per climb
function _climb_up(A::Matrix{Int}, config::ClimbConfig, ws::Workspace{T, F}, current_det::Real, best_det::Real; _cancel_flag) where {T<:Integer, F<:AbstractFloat}
    (; gradient, A_precise, A_inv, I_mat, tabu_mat, tabu_order_queue) = ws
    (; n_size, c_climbs, tabu_max, seed, verbose) = config

    if _cancel_flag[]
        return false, 0.0
    end

    move, delta, move_found = _best_change(A, tabu_mat, gradient, current_det, best_det)

    if move_found
        i, j = move.row, move.col
        A[i, j] = 1 - A[i, j]
        
        # maintain set of moves in bool matrix
        tabu_mat[i, j] = true
        # maintain order of moves
        push!(tabu_order_queue, (i, j)) 

        # grab oldest move and discard
        if length(tabu_order_queue) > tabu_max
            old_i, old_j = popfirst!(tabu_order_queue)
            tabu_mat[old_i, old_j] = false # Free the oldest Tabu move
        end

        return true, delta
    else
        return false, 0.0
    end
end


# Main climbing loop to perform many trials
# Called once to initiate sampling/simulation
function _hill_climb(config::ClimbConfig, t_trials::Int, ws::Workspace{T, F}; _cancel_flag) where {T<:Integer, F<:AbstractFloat}
    # Structural Setup
    # unpack config and workspace
    (; gradient, A_precise, A_inv, I_mat, tabu_mat, tabu_order_queue) = ws
    (; n_size, c_climbs, tabu_max, seed, verbose) = config

    # Other constructions for memory
    best_matrix = zeros(Int, n_size, n_size)
    prog_check = max(1, div(t_trials, 5)) # set progress check to every 20%
    shakes = Int(cld(n_size, 2))
    
    # Trackers
    best_det = typemin(T) # tabu metric
    max_climb = 0 # benchmarking
    current_trial = 0 # verbosity
    last_improvements = zeros(Int, t_trials) # benchmarking

    try
        for t in 1:t_trials

            
            current_trial = t
            # reset tabu data per trial
            fill!(tabu_mat, false)
            empty!(tabu_order_queue)

            # Allows for reactive interruptions
            if mod(t, 1) == 0
                yield() # Gives time for input (keyboard) check (negligible time)
            end
            if _cancel_flag[]
                break
            end

            # Reset starting matrix
            if isnothing(seed) || size(seed, 1) != n_size
                HCM = rand(0:1, n_size, n_size)
            else
                HCM = copy(seed)
                shuffle_bits(HCM, n_size, shakes, tabu_mat, tabu_order_queue, tabu_max)
            end

            # Find starting determinant
            if n_size < 23
                current_det = T(round(det(HCM))) # normal det() fine for smaller matrices
            else
                current_det = LinearAlgebra.det_bareiss(T.(HCM)) # det_bareiss technique needed for exactness in larger matrices
            end

            # Check if start > best
            if current_det > best_det
                best_det = current_det
                best_matrix = copy(HCM)
            end

            # Climbing Loop
            count = 0
            for c in 1:c_climbs
                # gradient updated in place below
                current_det = matrix_gradient(HCM, current_det, ws, tabu_max)
                moved, delta = _climb_up(HCM, config, ws, current_det, best_det; _cancel_flag=_cancel_flag)

                if !moved
                    count = c
                    break
                end

                # Update best
                current_det += T(delta)
                if current_det > best_det 
                    best_det = current_det
                    best_matrix = copy(HCM)
                end

                count = c 
            end 

            # Max climb tracker (for testing mostly)
            if count > max_climb
                max_climb = count
            end

            # Progress check
            if verbose
                if mod(t, prog_check) == 0
                    println("Core $(Threads.threadid()) completed Trial $t. Current Max Det: $(round(Int, abs(det(best_matrix))))")
                end
            end

        end

    catch err
        if err isa InterruptException
            println("\nEscape (Ctrl+C) detected. Exiting simulation early.")
            println("Core $(Threads.threadid()) exiting early at Trial $current_trial. Returning results...\n")
            _cancel_flag[] = true
        else
            rethrow(err)
        end
    end

    if verbose
	    println("Core $(Threads.threadid()) Simulation ended. Highest Climb Iteration: $max_climb. Trials completed: $current_trial")
    end

    return best_matrix
end


function _io_output(config, best_mat, final_time)
    (; n_size, c_climbs, t_trials, tabu_max, seed, verbose, parallel, thread_count, output) = config
    timestamp = round(BigInt, time())
    mkpath(joinpath(@__DIR__, "..", "results"))
    filepath = joinpath(@__DIR__, "..", "results", "C-$timestamp.txt")
    open(filepath, "w") do io
        println(io, "----------------------------------------------------------------")
        println(io ,"Timestamp: ", cfmt("%d", timestamp))
        println(io, "Size: $n_size x $n_size")
        println(io, "Determinant Value: ", cfmt("%'d", det(best_mat)))
        println(io, " --- Metrics & Parameters --- ")
        println(io, "Trials: $t_trials")
        println(io, "Max Climbs: $c_climbs")
        println(io, "Tabu Memory: $tabu_max")
        println(io, "Time taken: ", cfmt("%.2f", final_time), " seconds")
        println(io, "Parallel: $parallel")
        seeded = !isnothing(seed)
        println(io, "Seeded: $seeded")
        println(io, "----------------------------------------------------------------")
        println(io, "Best Matrix Found:")   
        show(IOContext(io, :limit => false), "text/plain", best_mat)
        println(io)
        println(io, "----------------------------------------------------------------")
    end
end


# Primary wrapper function for simple / intuitive use
function simulate_climb(config::ClimbConfig)
    (; n_size, c_climbs, t_trials, tabu_max, seed, verbose, parallel, thread_count, output) = config
    
    # For memory purposes. 
    # Big matrices also cause overflow/rounding concerns.
    if n_size < 23
        T, F = Int64, Float64
    elseif n_size <= 65
        T, F = Int128, Float64
    else
        T, F = BigInt, BigFloat
    end

    # Multithreaded
    if parallel
        LinearAlgebra.BLAS.set_num_threads(1)
        # Thread setup
        if thread_count < 1 || thread_count > Threads.nthreads() # safeguard for parameter
            thread_count = Threads.nthreads()
            println("Beginning multithreaded search using $(thread_count) threads. (Machine Detected)")
        else
            println("Beginning multithreaded search using $(thread_count) threads. (User Input)")
        end

        workspaces = [initialize_workspace(n_size, tabu_max, T, F) for _ in 1:thread_count]
        chunk_results = Vector{Matrix{Int64}}(undef, thread_count) # thread results location
        for i in 1:thread_count
            chunk_results[i] = zeros(Int, n_size, n_size)
        end
        c_flag = Threads.Atomic{Bool}(false) # interruption / early escape flag

        println("Total Trials: $t_trials (Distributed). Max Climbs Per Trial: $c_climbs. Max Tabu Memory: $tabu_max.")
        println("--------------------------------------------------------------------------------------")
        start_time = time()        
    
        # Each thread runs its portion of distributed trial load.
        try
            @threads for thread in 1:thread_count
                thread_space = workspaces[thread]
                thread_id = Threads.threadid()
                println("Core $thread_id starting new search on a $n_size x $n_size Matrix.")
                thread_trials = div(t_trials, thread_count)
                chunk_results[thread] = _hill_climb(config, thread_trials, thread_space; _cancel_flag=c_flag)
            end

        catch err
            if err isa InterruptException
                println("\nEscape (Ctrl+C) detected. Exiting simulation early.")
                println("Thread orchestrator winding down...")
                c_flag[] = true # signal to break
                sleep(0.5) # allow time for thread wind down
            else
                rethrow(err)
            end
        end

        # Find best matrix from all threads
        best_id = argmax(k -> det(chunk_results[k]), keys(chunk_results))
        best_mat = chunk_results[best_id]

    # Single threaded
    else
        c_flag = Threads.Atomic{Bool}(false)       

        println("Starting new search on a $n_size x $n_size Matrix.")
        println("Trials: $t_trials. Max Climbs Per Trial: $c_climbs. Max Tabu Memory: $tabu_max.")
        println("----------------------------------------------------------------")
        start_time = time()
        workspace = initialize_workspace(n_size, tabu_max, T, F)
        best_mat = _hill_climb(config, t_trials, workspace; _cancel_flag=c_flag)

    end

    # Display results from above process
    println("Best Matrix Found:")   
    show(IOContext(stdout, :limit => false), "text/plain", best_mat)
    println()

    best_det = round(BigInt, det(best_mat))
    println("Determinant Value: ", cfmt("%'d", best_det))
    final_time = time() - start_time
    println("Time taken: ", cfmt("%.2f", final_time), " seconds")

    if output
        # safe output file (incase terminal breaks)
        _io_output(config, best_mat, final_time)
    end
end

end # Module End

# hello reader!