println("allocating cores")
using Distributed, ClusterManagers
# automatically decide if it needs to be run in slurm envirnoment
try
    println("allocating $(ENV["SLURM_NTASKS"]) slurm tasks")
    addprocs_slurm(parse(Int64, ENV["SLURM_NTASKS"]))
catch err
    if isa(err, KeyError)
        N_tasts_ = 20
        println("allocating $N_tasts_ normal tasks")
        addprocs(N_tasts_)
    end
end

println("done")
flush(stdout); flush(stderr)

println("loading packages")
@everywhere using GadgetIO, GadgetUnits
@everywhere using Printf
@everywhere using Base.Threads
@everywhere using Statistics
@everywhere using ProgressMeter
@everywhere using DiffusiveShockAccelerationModels
println("done")
flush(stdout); flush(stderr)

# mapping settings
@everywhere const global snap_base = "/e/ocean3/Local/3072/nonrad_mhd_crs_new/snapdir_000_z=0/snap_000"

@everywhere global const GU = GadgetPhysical(GadgetIO.read_header(snap_base))

const jobs = RemoteChannel(()->Channel{Int}(32))
const results = RemoteChannel(()->Channel{Tuple}(1600))

# phase map settings
@everywhere const x_lim = [0.0, 3.0]
@everywhere const Nbins = 30

@everywhere include("bin_1D.jl")


@everywhere function dissipated_energy_1D_bin_of_subfile(subfile)

    println("subfile $subfile running on $(nthreads()) threads")
    flush(stdout)
    flush(stderr)

    data = Dict(block => read_block(snap_base * ".$subfile", block, parttype=0)
                for block ∈ ["MACH", "SHRH", "SHSP", "SHOB"])

    sel = findall(data["MACH"] .> 1.0)
    println("\tnumber of particles: $(length(sel))")

    mach = log10.(data["MACH"][sel])
    Fkin = @. 0.5 * data["SHRH"][sel] * data["SHSP"][sel]^3 * GU.rho_cgs * GU.v_cgs^3
    weights = ηB_acc_p.(data["SHOB"][sel])

    weights[isnan.(weights)] .= 0.0

    mach_count, Fkin_sum = bin_1D(mach, x_lim, Fkin, show_progress=false; Nbins)
    mach_count, Fkin_weighted_sum = bin_1D(mach, x_lim, weights .* Fkin, show_progress=false; Nbins)
    data = nothing
    GC.gc()

    return mach_count, Fkin_sum, Fkin_weighted_sum
end



@everywhere function write_Fkin_binning(filename, sum_mach_count, sum_Fkin, sum_Fkin_weighted)
    f = open(filename, "w")
    write(f, Nbins)
    write(f, x_lim)
    write(f, sum_mach_count)
    write(f, sum_Fkin)
    write(f, sum_Fkin_weighted)
    close(f)
end




@everywhere function do_work(jobs, results) # define work function everywhere
    while true
        job_id = take!(jobs)
        put!(results, dissipated_energy_1D_bin_of_subfile(job_id - 1))
    end
end

function make_jobs(n)
    for i in 1:n
        put!(jobs, i)
    end
end


function run_Fkin_1D()

    println("starting workers")

    for p in workers() # start tasks on the workers to process requests in parallel
        remote_do(do_work, p, jobs, results)
    end
    
    n = 2048
    #n = 4

    sum_mach_count    = zeros(Int64, Nbins)
    sum_Fkin          = zeros(Float64, Nbins)
    sum_Fkin_weighted = zeros(Float64, Nbins)


    errormonitor(@async make_jobs(n)); # feed the jobs channel with "n" jobs
    
    println("running")
    flush(stdout); flush(stderr)

    @time while n > 0 # print out results

        mach, Fkin, Fkin_weighted = take!(results)
        sum_mach_count    .+= mach
        sum_Fkin          .+= Fkin
        sum_Fkin_weighted .+= Fkin_weighted
        n -= 1

        mach = Fkin = Fkin_weighted = nothing
        GC.gc()
    end
    flush(stdout); flush(stderr)

    println("maximum = $(maximum(sum_Fkin))")
    flush(stdout); flush(stderr)
    filename = data_path * "Fkin.dat"
    write_Fkin_binning(filename, sum_mach_count, sum_Fkin, sum_Fkin_weighted)


end

run_Fkin_1D()