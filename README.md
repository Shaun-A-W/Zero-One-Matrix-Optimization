# Zero-One-Matrix-Optimization
Julia implementation for Hill Climbing and Targeted Search algorithms on Zero One Matrix Determinants

## Dependencies

Refer to .toml files for dependencies/requirements. Utilize Julia PKG for ease of handling this. Some methods below:

### Inside Julia REPL

1) Start Julia inside project directory. 
2) Enter pkg mode by pressing "]".
3) Activate & Download via:
    
    pkg> activate .
    
    pkg> instantiate

### From Command Line/Terminal

Navigate to directory and enter the following:

julia --project=. -e 'using Pkg; Pkg.instantiate()'


## Direction

Made to be readable, though some optimization efforts to reduce the strain and use of the Julia Garbage Collector may make some functions confusing. Start with docstrings and comments in the flow of function calls for best intuition.

Refer to /src/ folder for further documentation on functionality and structure. /src/README.md contains summarized instruction, but I endorse spending time understanding the program's workflow.