using Plots
using LinearAlgebra
using Format
using Statistics
using StatsBase

function new_matrix(n_size::Int)
    return rand(0:1, n_size, n_size)
end


function _data_increment!(change::Tuple{Int, Int}, data::Dict{Tuple{Int, Int}, Int})
    data[change] = get!(data, change, 0) + 1
end


function _markov_collection(n_size::Int, n_seeds::Int, n_flips::Int)
    if n_seeds < 1 || n_flips < 1
        throw(ValueError)
    end

    # store changes as (old, new) , count
    _data_dump = Dict{Tuple{Int, Int}, Int}()

    n_elem = n_size^2

    for origin in 1:n_seeds
        seed = new_matrix(n_size)
        seed_float = Float64.(seed)
        det_old = det(seed_float)

        for flip in 1:n_flips
            # perform random flip
            r_coord = rand(1:n_elem)
            seed[r_coord] = xor(seed[r_coord], 1) # faster since binary 0s & 1s
            seed_float[r_coord] = 1.0 - seed_float[r_coord]

            # record change
            det_new = det(seed_float)
            change = (round(Int, det_old), round(Int, det_new))
            _data_increment!(change, _data_dump)

            det_old = det_new
        end
    end
    
    return _data_dump
end


function _markov_scatter(n_size, befores, afters, counts, max_count, color, caption)
    # normalize and scale sizes
    scaled_sizes = 3 .+ (counts ./ max_count) .* 30 # base + proportion * max - base

    # scale opacity dynamically
    scaled_alpha = 0.05 .+ 0.10 .* (counts ./ max_count)

    markov_plot = scatter(
        befores, 
        afters, 
        markersize = scaled_sizes, 
        markeralpha = scaled_alpha,
        color = color,
        markerstrokewidth = 0.8, # subtle outlines help distinguish overlapping bubbles
        legend = false,
        title = "\nDeterminant Transition States ($n_size x $n_size)",
        xlabel = "\nPre-flip Determinant (Old)\n$caption",
        ylabel = "Post-flip Determinant (New)\n",
        grid = true,
        size = (700,600),
        top_margin = 5Plots.mm,
        bottom_margin = 5Plots.mm, 
        left_margin = 5Plots.mm,
        right_margin = 10Plots.mm
    )

    display(markov_plot)
end


function _markov_heat(n_size, befores, afters, counts, theme, caption)
    # # find how many unique values exist along each axis
    # uniq_x_cnt = length(unique(befores))
    # uniq_y_cnt = length(unique(afters))

    # # dynamic target: sqrt of unique counts, clamped between values
    # x_bins = clamp(round(Int, sqrt(uniq_x_cnt)), 100, 1000)
    # y_bins = clamp(round(Int, sqrt(uniq_y_cnt)), 100, 1000)

    # # for a perfectly square aspect ratio grid, pick the max of the two:
    # # dense option  (grouped)
    # final_nbins = max(x_bins, y_bins)
    # println(final_nbins)

    # # sparse option (granular)
    # # total_unique_states = length(unique([befores; afters]))
    
    # Freedman-Diaconis Method
    # Compute the robust spread (IQR) of your data
    iqr_x = quantile(befores, 0.75) - quantile(befores, 0.25)
    width_x = 2 * iqr_x * length(befores)^(-1/3)

    # Calculate exactly how many bins fit into your data span
    span_x = maximum(befores) - minimum(befores)
    final_nbins = round(Int, span_x / width_x)

    h2d = fit(Histogram, (befores, afters), weights(counts), nbins=max(final_nbins, n_size*5))

    edges_x = h2d.edges[1]  # Binned boundaries for Pre-flip
    edges_y = h2d.edges[2]  # Binned boundaries for Post-flip
    raw_grid = h2d.weights  # The aggregated 2D count matrix

    log_matrix = [c > 0 ? log10(c) : NaN for c in raw_grid]

    heat_map_log = heatmap(
        edges_x,
        edges_y,
        log_matrix,
        cmap = theme,
        clims = (0, maximum(filter(!isnan, log_matrix))),
        title = "\n(Log10) Determinant Transition States ($n_size x $n_size)",
        colorbar_title = "Log10(Counts)",
        xlabel = "\nPre-flip Determinant (Old)\n$caption",
        ylabel = "Post-flip Determinant (New)\n",
        grid = true,
        size = (700,600),
        top_margin = 5Plots.mm,
        bottom_margin = 5Plots.mm, 
        left_margin = 5Plots.mm,
        right_margin = 10Plots.mm
    )

    display(heat_map_log)
end


function markov_routine(n_size::Int, n_seeds::Int, n_flips::Int, plot_type::String; scatter_color::Symbol=:purple, heat_theme::Symbol=:turbo)
    @assert n_size > 0 "Parameter Error: n_size must be greater than 0 (got $n_size)"
    @assert n_seeds > 0 "Parameter Error: n_seeds must be greater than 0 (got $n_seeds)"
    @assert n_flips > 0 "Parameter Error: n_flips must be greater than 0 (got $n_flips)"

    if plot_type != "scatter" && plot_type != "heat"
        error("Parameter Error: Invalid plot type ($plot_type). Supported types : 'scatter' , 'heat'")
    end

    println("\nStarting Markov Chain Simulation on $n_size dimensional matrices...")
    start = time()
    results = _markov_collection(n_size, n_seeds, n_flips)

    changes = collect(keys(results))
    counts = collect(values(results))
    max_count = maximum(counts)

    n_unique = length(results)
    eff_flips = sum(counts)
    caption_text = "Number of Unique Changes: $n_unique \nNumber of Effective Flips: $(format(eff_flips, commas=true))"

    befores = [c[1] for c in changes] # pre-flip determinants
    afters = [c[2] for c in changes] # post-flip determinants

    if plot_type == "scatter"
        _markov_scatter(n_size, befores, afters, counts, max_count, scatter_color, caption_text)
    end

    if plot_type == "heat"
        _markov_heat(n_size, befores, afters, counts, heat_theme, caption_text)
    end

    println("Finished analysis. See plot for relative frequency of changes.")
    println("Time Taken: ", cfmt("%.2f", (time() - start)), "s")
end


# note the following max determinant values
# max_determinants = Dict(
#     1 => 1,
#     2 => 1,
#     3 => 2,
#     4 => 3,
#     5 => 5,
#     6 => 9,
#     7 => 32,
#     8 => 56,
#     9 => 144,
#     10 => 320,
#     11 => 1458,
#     12 => 3645,
#     13 => 9477,
#     14 => 25515,
#     15 => 131072,
#     16 => 327680,
#     17 => 1114112,
#     18 => 3411968,
#     19 => 19531250,
#     20 => 56640625,
#     21 => 195312500,
#     22 => 967396000,
# )

# nice themes:
# :turbo
# :heat
# :blues
# cgrad(:magma, rev=true)

# n_size = 12
# n_seeds = 10_00
# n_flips = 10_00
# plot_type = "heat"

# markov_routine(n_size, n_seeds, n_flips, plot_type)

for n in 5:15
    markov_routine(n, 10_0, 10_00000, "heat")
end