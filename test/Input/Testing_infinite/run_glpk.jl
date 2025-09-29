## This file should be in a directory containing an Input directory, which has a directory
# called <data_dir> containing the JADE input files.

using JADE, JuMP

# Import the GLPK solver
using GLPK 
optimizer = GLPK.Optimizer  


## Set directory containing the Input / Output subdirectories
cd(@__DIR__)
ENV["JADE_DIR"] = dirname(dirname(@__DIR__))
@eval JADE JADE_DIR = dirname(dirname(@__DIR__))

## Modify these settings
data_dir = basename(@__DIR__)
runfile = "run"
training = true
simulation = false

## Set up inputs for JADE models
rundata = define_JADE_model(data_dir, run_file = runfile)

# ## Create JADE model from the runfile data
model = create_JADE_model(rundata, optimizer)

if training
    solve_options = define_JADE_solve_options(data_dir, run_file = runfile)
    optimize_policy!(model, solve_options)
end

if simulation
    sim_settings = define_JADE_simulation(data_dir, run_file = runfile)
    results = simulate(model, sim_settings)
end
