Direct self to program files to read function comments and docstrings for explanations and examples.

## Usage

It is recommended to utilize both programs in separate files, like the main.jl example files listed, to avoid conflicts or issues. As found in docstrings and example main_climb.jl and main_target.jl testing files, the structure of the user interface is simple. Calling these two interfaces looks like this:

### Hill Climbing Algorithm

Configuration Example

```julia
include("HillClimbZeroOne.jl")
using .HillClimbZeroOne

my_config = ClimbConfig(
	n_size = 22,			# n x n matrix.
	c_climbs = 30000,		# Max climbs per trial.
	t_trials = 100,		    # Total trials (random starts). Distributed if parallel.
	tabu_max = 30,			# Max memory for tabu list.
	parallel = false,       # Boolean to determine thread mode.
	verbose = true,         # Boolean to determine extent of console output.
	seed = nothing,         # Can take a matrix as a seed (shuffled slightly) instead of random per trial.
	thread_count = 0,       # Determines the number of threads to use. Enter 0 to automatically determine.
	output = true           # Boolean to determine whether to output a results file.
)

HillClimbZeroOne.simulate_climb(my_config)
```

### Target Search Algorithm

Configuration Example

```julia
include("TargetClimbZeroOne.jl")
using .TargetClimbZeroOne

my_config = TargetConfig(
	target = 119993         # Target determinant value.
	n_size = 22,			# n x n matrix.
	c_climbs = 30000,		# Max climbs per trial.
	t_trials = 100,		    # Total trials (random starts). Distributed if parallel.
	tabu_max = 30,			# Max memory for tabu list.
	parallel = false,       # Boolean to determine thread mode.
	verbose = true,         # Boolean to determine extent of console output.
	seed = nothing,         # Can take a matrix as a seed (shuffled slightly) instead of random per trial.
	thread_count = 0,       # Determines the number of threads to use. Enter 0 to automatically determine.
	output = true           # Boolean to determine whether to output a results file.
)

TargetClimbZeroOne.simulate_target(my_config)
```