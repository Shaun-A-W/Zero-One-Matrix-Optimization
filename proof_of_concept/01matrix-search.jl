# Shaun W.
# UTulsa Department of Mathematics
# Analysis of Zero-One Matrices
# Optimization / Search Algorithms


using LinearAlgebra
using Format
using Base.Threads


function matrix_gradient(A::Matrix{Int})
    n = size(A, 1)
	d = det(A)
	
	# If singular, use general method
    if abs(d) < 1e-9
        grad = zeros(Int, n, n)
        for i in 1:n
            for j in 1:n
                # Using @views avoids memory allocations during slicing
                sub_matrix = @views A[1:n .!= i, 1:n .!= j]
                grad[i, j] = ((-1)^(i + j)) * round(Int, det(sub_matrix))
            end
        end
        return grad
    end
    
    # If invertible, use shortcut
    return round.(Int, d .* inv(A)')
end


function best_change(A::Matrix{Int}, tabu_list::Vector{Tuple{Int, Int}})
	n = size(A, 1)
	gradient = matrix_gradient(A)
    fallback_move = (row = 0, col = 0)
    fallback_delta = -Inf
    fallback_found = false

	# Note: Flipping 0 -> 1 adds grad
	#		Flipping 1 -> 0 subs grad

	 for i in 1:n 
        for j in 1:n
            is_tabu = (i, j) in tabu_list
            delta = A[i, j] == 0 ? gradient[i, j] : -gradient[i, j]
            
            if !is_tabu
                # First-Improvement: If we find a strictly positive move
                if delta > 0
                    return (row = i, col = j), delta, true
                end
                
                # Else, track the least-bad non-tabu move to navigate flat plateaus
                if delta > fallback_delta
                    fallback_delta = delta
                    fallback_move = (row = i, col = j)
                    fallback_found = true
                end
            end
        end
    end
	return fallback_move, fallback_delta, fallback_found # best improvement / last resort 
end


function target_change(A::Matrix{Int}, tabu_list::Vector{Tuple{Int, Int}}, target::Int)
	# Unlike best_change(),
	# Try to hit specific determinant target
	n = size(A, 1)
	d = det(A)
	gradient = matrix_gradient(A)
    fallback_move = (row = 0, col = 0)
    fallback_proximity = abs(target - abs(d))
    fallback_found = false

	# Note: Flipping 0 -> 1 adds grad
	#		Flipping 1 -> 0 subs grad

	 for i in 1:n 
        for j in 1:n
            is_tabu = (i, j) in tabu_list
            delta = A[i, j] == 0 ? gradient[i, j] : -gradient[i, j]
			new_d = d + delta
			proximity = abs(target - abs(new_d))
            
            if !is_tabu
				if proximity <= fallback_proximity
					fallback_proximity = proximity
					fallback_move = (row = i, col = j)
					fallback_found = true
				end
            end
        end
    end
	return fallback_move, fallback_proximity, fallback_found # best improvement / last resort 
end


function do_best(A::Matrix{Int}, tabu_list::Vector{Tuple{Int, Int}}, tabu_max, target=nothing)
	if isnothing(target)
		move, change, move_found = best_change(A, tabu_list)

		if !move_found
			return false
		end

		i, j = move.row, move.col
		A[i, j] = 1 - A[i, j] # bit flip

		# Update Tabu
		push!(tabu_list, (i, j))
		if length(tabu_list) > tabu_max
			popfirst!(tabu_list) # Remove oldest
		end

		return true

	else
		move, proximity, move_found = target_change(A, tabu_list, target)

		if move_found

			i, j = move.row, move.col
			A[i, j] = 1 - A[i, j] # bit flip

			# Update Tabu
			push!(tabu_list, (i, j))
			if length(tabu_list) > tabu_max
				popfirst!(tabu_list) # Remove oldest
			end
			return true, proximity

		else
			return false, Inf
		end
	end
end


function main_loop(n_size, c_climbs, t_trials, tabu_max, tabu_list)
	best_matrix = zeros(Int, n_size, n_size)
	max_climb = 0
	first_run = true

	for t in 1:t_trials
		# Perform hill climb on random matrix
		HCM = rand([0,1], n_size, n_size)
		tabu_list = Tuple{Int, Int}[] # reset tabu memory

		count = 0
		for c in 1:c_climbs
			if !do_best(HCM, tabu_list, tabu_max)
				# println("Stopping loop at move $c")
				# println("Gradient of Changed A:")
				# display(matrix_gradient(HCM))
				
				count = c
				break
			end
			count = c
		end

		if count > max_climb
			max_climb = count
		end

		# Check for new best
		if first_run || abs(det(HCM)) > abs(det(best_matrix))
			best_matrix = copy(HCM)
			first_run = false
		end

		# Progress check
		if mod(t, 500) == 0
			println("Trial $t Completed on Core $(threadid()). Current Thread Max Det: $(round(Int, abs(det(best_matrix))))")
		end

	end
	println("Highest Climb Loop: $max_climb")
	return best_matrix
end


function main_loop_target(n_size, c_climbs, t_trials, tabu_max, tabu_list, target)
    best_matrix = zeros(Int, n_size, n_size)
    best_proximity = Inf
    target_hit = false
	end_move = nothing

    for t in 1:t_trials
        HCM = rand(0:1, n_size, n_size)
        tabu_list = Tuple{Int, Int}[] 
        
        # Track the best proximity in this specific trial
        trial_best_prox = abs(target - abs(det(HCM)))

        for c in 1:c_climbs
            success, proximity = do_best(HCM, tabu_list, tabu_max, target)
            
            if !success
                break
            end

            if proximity < trial_best_prox
                trial_best_prox = proximity
                # println("Trial $t | Move $c | Distance to Target dropped to: $(cfmt("%'d", proximity))")
            end

            if proximity < 0.5
                # println("TARGET FOUND ON TRIAL $t, MOVE $c")
                best_matrix = copy(HCM)
                target_hit = true
				end_move = c
                break
            end
        end

        # Break out if target found
        if target_hit
			best_proximity = 0
			best_matrix = copy(HCM)
            break
        end

		# Check for new best
        if trial_best_prox < best_proximity
            best_proximity = trial_best_prox
            best_matrix = copy(HCM)
        end

		# Progress check
        if mod(t, 100) == 0
            prox_int = isinf(best_proximity) ? 0 : round(Int, best_proximity)
            if prox_int == 0 && !isinf(best_proximity)
                # If it rounds to 0 but didn't trigger a break, it's floating-point noise
                prox_string = "Near Match (<1)"
            else
                prox_string = isinf(best_proximity) ? "Infinity" : cfmt("%'d", prox_int)
            end
            # println("Trial $t Completed. Closest proximity so far: $prox_string")
        end
    end
	if best_proximity != 0
		# println("TARGET NOT FOUND")
	end
    return best_matrix, best_proximity, end_move
end


# Main block to avoid unintended executions
if !any(frame -> frame.func == :include, stacktrace())

	####################
	# Variables & Main #
	####################

	n_size = 21				# n x n matrix
	c_climbs = 1_000		# max climbs per trials
	t_trials = 1_000		# total trials (random starts)
	tabu_max = 10			# max memory for tabu list
	tabu_list = Tuple{Int, Int}[]	# empty tabu list, do not touch
	target = 195312500
	# for parallel only...
	chunks = 20
    loops_per_chunk = 1_000


	# HILL CLIMBING SEARCH
	# start_time = time()
	# println("Starting new search on a $n_size x $n_size Matrix.")
	# best = main_loop(n_size, c_climbs, t_trials, tabu_max, tabu_list)

	# println("Best Matrix Found:")
	# show(IOContext(stdout, :limit => false), "text/plain", best)
	# println()

	# best_det = round(Int, det(best))
	# println("Determinant Value: ", cfmt("%'d", best_det))

	# end_time = time()
	# println("Time taken: ", cfmt("%.2f", end_time - start_time), " seconds")



	# TARGET SEARCH -- INCLUDES NEGATIVE PAIR (I.E. MAY RETURN NEGATIVE TARGET ... ROW SWAPS)
	# start_time = time()
	# println("Starting new search for $target on a $n_size x $n_size Matrix.")
	# best_target, best_prox, end_move = main_loop_target(n_size, c_climbs, t_trials, tabu_max, tabu_list, target)

	# println("Closest Matrix Found:")
	# show(IOContext(stdout, :limit => false), "text/plain", best_target)
	# println()

	# best_det = round(Int, det(best_target))
	# println("Determinant Value: ", cfmt("%'d", best_det))
	# println("Proximity: ", cfmt("%'d", best_prox))

	# end_time = time()
	# println("Time taken: ", cfmt("%.2f", end_time - start_time), " seconds")



	# PARALLEL HILL-CLIMBING
    # Array to hold the resulting matrices
    chunk_results = Vector{Matrix{Int}}(undef, chunks)
    
    # Run the chunks in parallel
	start_time = time()
    @threads for c in 1:chunks
        # Each thread runs main_loop for 10,000 internal trials
		println("Core $(threadid()) starting new search on a $n_size x $n_size Matrix.")
		tabu_mem = (threadid() * 2)
        chunk_results[c] = main_loop(n_size, c_climbs, loops_per_chunk, tabu_mem, tabu_list)
    end
	
	best_id = argmax(k -> det(chunk_results[k]), keys(chunk_results))
	println("Best matrix from multithreaded search:")
	show(IOContext(stdout, :limit => false), "text/plain", chunk_results[best_id])
	println()
	println("Determinant Value: ", cfmt("%'d", det(chunk_results[best_id])))
	end_time = time()
	println("Time taken: ", cfmt("%.2f", end_time - start_time), " seconds")

end