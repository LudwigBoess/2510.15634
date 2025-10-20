println("allocating cores")
using Distributed, ClusterManagers

# automatically decide if it needs to be run in slurm envirnoment
try
    println("allocating $(ENV["SLURM_NTASKS"]) slurm tasks, using $(2 * parse(Int64, ENV["SLURM_CPUS_PER_TASK"])) threads each")
    withenv("JULIA_NUM_THREADS" => 2 * parse(Int64, ENV["SLURM_CPUS_PER_TASK"])) do
        addprocs_slurm(parse(Int64, ENV["SLURM_NTASKS"])) # spawn 4 workers with 2 threads each
    end
catch err
    if isa(err, KeyError)
        println("allocating 4 normal tasks")
        addprocs(4)
    end
end

println("loading packages")
@everywhere using GadgetIO, GadgetUnits
@everywhere using Printf
@everywhere using SPHKernels, SPHtoGrid
@everywhere using SpectralCRsUtility
@everywhere using Base.Threads
@everywhere using Statistics
@everywhere using Healpix

println("done")
flush(stdout);
flush(stderr);


# include map functions
@everywhere include("config.jl")
@everywhere include("crpp.jl")
@everywhere include("gamma.jl")

"""
    run_it()

Main function
"""
function run_it()
    if ARGS[1] == "CRpP"
        filename = map_path * "allsky_CRpP_$viewpoint.fits"
        # `distributed_allsky_map` is imported from SPHtoGrid.jl
        distributed_allsky_map(filename, Nside, Nfiles, CRpP_maps_of_subfile)
    elseif ARGS[1] == "Xcr"
        filename = map_path * "allsky_Xcr_$viewpoint.fits"
        distributed_allsky_map(filename, Nside, Nfiles, Xcr_maps_of_subfile)
    elseif ARGS[1] == "Xcr_mean"
        filename = map_path * "allsky_Xcr_mean_$viewpoint.fits"
        distributed_allsky_map(filename, Nside, Nfiles, Xcr_mean_maps_of_subfile)
    elseif ARGS[1] == "gamma"
        filename = map_path * "allsky_gammaSB_$viewpoint.fits"
        distributed_allsky_map(filename, Nside, Nfiles, gamma_maps_of_subfile, reduce_image=false)
    end
end

run_it()