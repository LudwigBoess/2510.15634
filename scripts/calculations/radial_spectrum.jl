using GadgetIO, GadgetUnits
using SpectralCRsUtility
using Base.Threads
using Printf
using ProgressMeter
#using LoopVectorization
#using BenchmarkTools
using PyPlot, PyPlotUtility
using DelimitedFiles
using PyCall

const global snap_base = "/e/ocean3/Local/3072/nonrad_mhd_crs_new/snapdir_000_z=0/snap_000"
const global sub_base = "/e/ocean3/Local/3072/nonrad_mhd_crs_new/groups_000_z=0/sub_000"

function get_r(pos, gpos)

    Npart = length(pos[1, :])
    r = Vector{Float64}(undef, Npart)

    @inbounds for i = 1:Npart
        r[i] = 0.0
        for dim = 1:3
            r[i] += (pos[dim, i] - gpos[dim])^2
        end
        r[i] = √(r[i])
    end

    return r
end


"""
    get_CRpE(data)

Energy spectrum in GeV
"""
function get_CRpE(data)

    Nbins = size(data["CRpN"], 1)
    par = CRMomentumDistributionConfig(0.1, 1.e5, Nbins)
    bounds = momentum_bin_boundaries(par)

    GU = GadgetPhysical(read_header(snap_base))
    Npart = length(data["RHO"])

    CRpE = Matrix{Float64}(undef, Nbins, Npart)

    @showprogress for i = 1:Npart

        norm = 10.0 .^ data["CRpN"][:, i]
        slope = Float64.(data["CRpS"][:, i])
        cut = Float64(data["CRpC"][i])
        rho = data["RHO"][i] * GU.rho_physical

        CRpE[:, i] .= energy_spectrum(norm, slope, cut, rho, bounds)
    end

    return CRpE .* data["MASS"][1]
end

function bin_data(r, CRpE, Nbins)

    bin_lim = [-2.0, log10(2.0)]
    # get logarithmic bin spacing
    dbin = (bin_lim[2] - bin_lim[1]) / Nbins

    logr = log10.(r)

    CRpE_bins = zeros(6, Nbins)
    CRpE_count = zeros(Int64, 6, Nbins)

    @showprogress for i = 1:length(r)
        if isinf(logr[i])
            continue
        end
        bin = 1 + floor(Int64, (logr[i] - bin_lim[1]) / dbin)
        # check if bin is in range to check if particle is relevant
        if (1 <= bin <= Nbins)

            # loop over all CR bins
            for cr_bin = 1:6
                # sum up logarithmic contribution of bin
                if (CRpE[cr_bin, i] > 0.0) && !isinf(CRpE[cr_bin, i])
                    # sum up to total binned quantity
                    CRpE_bins[cr_bin, bin] += log10(CRpE[cr_bin, i])
                    # count up histogram storage
                    CRpE_count[cr_bin, bin] += 1
                end
            end
        end
    end

    # get log mean
    @threads for i ∈ eachindex(CRpE_count)
        if !iszero(CRpE_count[i])
            CRpE_bins[i] /= CRpE_count[i]
        end
    end

    r_bin_centers = [bin_lim[1] + (i - 0.5) * dbin for i = 1:Nbins]

    return r_bin_centers, CRpE_bins
end

function get_data(halo_id, Nbins)

    blocks = ["POS", "MASS", "RHO", "CRpN", "CRpS", "CRpC"]

    @info "Reading subfind"
    # only in r500
    gpos = read_subfind(sub_base * ".$(halo_id.file)", "GPOS")[:, halo_id.id]
    r200 = read_subfind(sub_base * ".$(halo_id.file)", "R200")[halo_id.id]

    @info "reading data"
    sphere = GadgetSphere(gpos, 2r200)
    data = read_particles_in_geometry(snap_base, blocks, sphere, verbose=true)

    @info "get r"
    r = get_r(data["POS"], gpos) ./ r200

    @info "get CRpE"
    CRpE = get_CRpE(data)

    @info "binning"
    bin_data(r, CRpE, Nbins)
end

# Coma
halo_id = HaloID(1, 9)
Nbins = 10

#r_bins, CRpE_bins = get_data(halo_id, Nbins)