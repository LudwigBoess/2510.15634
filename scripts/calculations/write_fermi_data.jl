"""
    Reproduce Fig. 05 from Ha+2020: 
    https://ui.adsabs.harvard.edu/abs/2020ApJ...892...86H/abstract

    Based on data by Ackermann+14:
    https://iopscience.iop.org/article/10.1088/0004-637X/787/1/18
    
    Gamma energy range: 500MeV - 200GeV
"""

using GadgetIO, GadgetUnits
using SpectralCRsUtility
using Statistics
using ProgressMeter
using Base.Threads
using Statistics
using DelimitedFiles
using Unitful, UnitfulAstro
using Cosmology
using Printf

struct ClusterGamma
    name::String
    M500::Float64
    Xcr::Float64
    Sγ::Float64
end

ClusterGamma(name, M::Float64, T::Tuple) = ClusterGamma(name, M, T[1], T[2])

function get_gamma_photons(data, h;
    pmin=0.1, pmax=1.e5)

    center = [247.980, 245.480, 255.290] .* 1.e3

    GU = GadgetPhysical(h)

    # cr setup 
    Nbins = size(data["CRpN"], 1)
    par = CRMomentumDistributionConfig(pmin, pmax, Nbins)
    bounds = [par.pmin * 10.0^((i - 1) * par.bin_width) for i = 1:par.Nbins+1]

    Npart = size(data["CRpN"], 2)

    Sγ = Vector{Float64}(undef, Npart)
    Lγ = Vector{Float64}(undef, Npart)

    @info "gamma flux"

    P = Progress(Npart)

    @threads for i = 1:Npart

        norm = GU.CR_norm .* 10.0 .^ data["CRpN"][:, i]
        slope = Float64.(data["CRpS"][:, i])
        cut = Float64(data["CRpC"][i])

        pos = data["POS"][:,i] - center 

        d = √(pos[1]^2 + pos[2]^2 + pos[3]^2) * GU.x_cgs

        nH = data["RHO"][i] * GU.rho_ncm3

        V = data["MASS"][i] / data["RHO"][i] * GU.x_cgs^3
        
        Sγ[i] = gamma_flux_pions(norm, slope, cut, bounds, nH, V, d, 
                                N_integration_steps=50,
                                Eγ_min=0.5, Eγ_max=200.0)

        Lγ[i] = Sγ[i] * 4π * d^2
        next!(P)
    end

    return sum(Sγ), sum(Lγ)
end

function get_data(halo_id::HaloID)

    sub_base = "/e/ocean3/Local/3072/nonrad_mhd_crs_new/groups_000_z=0/sub_000"
    snap_base = "/e/ocean3/Local/3072/nonrad_mhd_crs_new/snapdir_000_z=0/snap_000"

    blocks = ["POS", "MASS", "RHO", "U", "CRpP", "CRpN", "CRpS", "CRpC"]
    
    # only in r500
    gpos = read_subfind(sub_base * ".$(halo_id.file)", "GPOS")[:,halo_id.id]
    r200 = read_subfind(sub_base * ".$(halo_id.file)", "R200")[halo_id.id]

    sphere = GadgetSphere(gpos, r200)
    data = read_particles_in_geometry(snap_base, blocks, sphere, verbose=true)
    
    h = read_header(snap_base)
    GU = GadgetPhysical(h)

    @info "Xcr"
    γ_m1 = 5 / 3 - 1
    Pth = γ_m1 .* data["RHO"] .* data["U"] .* GU.P_th_cgs
    Pcr = data["CRpP"] .* GU.P_CR_cgs
    Xcr_sum = sum(Pcr) / sum(Pth)

    sel = findall(Pcr .> 0)
    Xcr_mean = mean(Pcr[sel] ./ Pth[sel])

    Sγ, Lγ = get_gamma_photons(data, h)

    return Xcr_sum, Xcr_mean, Sγ, Lγ
end


function get_M200(halo_id)

    sub_base = "/e/ocean3/Local/3072/nonrad_mhd_crs_new/groups_000_z=0/sub_000"
    h = read_header(sub_base)
    GU = GadgetPhysical(h)
    read_subfind(sub_base * ".$(halo_id.file)", "M200")[halo_id.id] * GU.m_msun |> Float64
end


function write_data(halo_id::HaloID, store_data_path::String)

    @info "M200"
    M200 = get_M200(halo_id)

    @info "reading data"
    Xcr_sum, Xcr_mean, Sγ, Lγ = get_data(halo_id)

    open(store_data_path * "gamma_R200.txt", "a") do io
        write(io, "Cluster_$(halo_id.file)_$(halo_id.id)\t$M200\t$Xcr_sum\t$Xcr_mean\t$Sγ\t$Lγ\n")
    end
end



# clusters = Vector{ClusterGamma}(undef, length(cluster_ids_Mgt5e4))
# @showprogress for i ∈ 1:length(cluster_ids_Mgt5e4)
#     cluster = cluster_ids_Mgt5e4[i]
#     clusters[i] = ClusterGamma("Cluster_$(cluster.file)_$(cluster.id)", get_M200(cluster), get_data(cluster))
# end

# println(clusters)


# const mp_c2 = 1.6726e-24^3 * 2.9979e10^4
# const GeVtoerg = 1.60217733e-3

function convert_ackermann(data_path)

    ackermann = readdlm(data_path * "ackermann.txt", ',', comments=true)
    m = ackermann[:, 2] .* 0.7 ./ 0.677

    c = cosmology(h=0.7, OmegaM=0.3)

    # gamma luminosity [1/s]
    Lγ = Vector{Float64}(undef, length(ackermann[:, 1]))

    for i = 1:length(ackermann[:, 1])
        z = ackermann[i, end]
        dA = angular_diameter_dist(c, z) * (1 + z) |> u"cm" |> ustrip
        Lγ[i] = ackermann[i, 5] * 4π * dA^2
    end

    return m, Lγ
end

# m, Lγ = convert_ackermann(data_path)


function get_data(snap_base::String)

    blocks = ["POS", "MASS", "RHO", "U", "CRpP", "CRpN", "CRpS", "CRpC"]
    
    # only in r500
    gpos = [245040.453125, 327781.84375, 246168.6875]
    r200 = 2642.8608

    sphere = GadgetSphere(gpos, r200)
    data = read_particles_in_geometry(snap_base, blocks, sphere, 
                                      use_keys=false, verbose=false)
    
    h = read_header(snap_base)
    GU = GadgetPhysical(h)

    @info "Xcr"
    γ_m1 = 5 / 3 - 1
    Pth = γ_m1 .* data["RHO"] .* data["U"] .* GU.P_th_cgs
    Pcr = data["CRpP"] .* GU.P_CR_cgs
    Xcr_sum = sum(Pcr) / sum(Pth)

    sel = findall(Pcr .> 0)
    Xcr_mean = mean(Pcr[sel] ./ Pth[sel])

    Sγ, Lγ = get_gamma_photons(data, h)
    Sγ_shifted, Lγ_shifted = get_gamma_photons(data, h, pmin=0.01, pmax=1.e4)

    println(Sγ, " ", Sγ_shifted)

    return Xcr_sum, Xcr_mean, Sγ, Sγ_shifted, Lγ, Lγ_shifted
end


function write_data(sim_name::String, snap_base::String, store_data_path::String)

    @info "M200"
    M200 = 1.945080326840785e15

    Xcr_sum, Xcr_mean, Sγ, Sγ_shifted, Lγ, Lγ_shifted = get_data(snap_base)

    open(store_data_path * "gamma_zoom_R200_selected.txt", "a") do io
        write(io, "$sim_name\t$M200\t$Xcr_sum\t$Xcr_mean\t$Sγ\t$Sγ_shifted\t$Lγ\t$Lγ_shifted\n")
    end
end

function convert_halo_ids(filename)

    data = readdlm(filename)
    cluster_ids = Vector{HaloID}(undef, size(data, 1))
    for i = 1:size(data,1)
        cluster_ids[i] = HaloID(data[i,1], data[i,2])
    end

    return cluster_ids
end


# cluster_ids_Mgt5e4 = [HaloID(0, 1), HaloID(0, 2), HaloID(0, 3), HaloID(0, 4), HaloID(0, 5), HaloID(0, 6),
#     HaloID(1, 1), HaloID(1, 2), HaloID(1, 3), HaloID(1, 4), HaloID(1, 5), HaloID(1, 6),
#     HaloID(1, 7), HaloID(1, 8), HaloID(1, 9), HaloID(1, 10), HaloID(1, 11), HaloID(1, 12),
#     HaloID(1, 13), HaloID(1, 14), HaloID(1, 15), HaloID(1, 16), HaloID(1, 17), HaloID(1, 18), HaloID(1, 19),
#     HaloID(2, 1), HaloID(2, 2), HaloID(2, 3), HaloID(2, 4), HaloID(2, 5), HaloID(2, 6), HaloID(2, 7),
#     HaloID(2, 8), HaloID(2, 9), HaloID(2, 10), HaloID(2, 11), HaloID(2, 12), HaloID(2, 13), HaloID(2, 14),
#     HaloID(2, 15), HaloID(2, 16), HaloID(2, 17), HaloID(2, 18), HaloID(2, 19), HaloID(2, 20), HaloID(2, 21),
#     HaloID(2, 22), HaloID(2, 23), HaloID(2, 24), HaloID(2, 25), HaloID(2, 26), HaloID(2, 27), HaloID(2, 28),
#     HaloID(2, 30), HaloID(2, 31), HaloID(2, 32), HaloID(2, 33), HaloID(3, 2), HaloID(3, 3), HaloID(3, 4),
#     HaloID(3, 5), HaloID(3, 7), HaloID(3, 8), HaloID(3, 9), HaloID(3, 10), HaloID(3, 11), HaloID(3, 12),
#     HaloID(3, 13), HaloID(3, 14), HaloID(3, 15), HaloID(3, 16), HaloID(3, 17), HaloID(3, 18), HaloID(3, 19),
#     HaloID(3, 20), HaloID(3, 21), HaloID(3, 22), HaloID(3, 25), HaloID(3, 27), HaloID(3, 29), HaloID(3, 30),
#     HaloID(3, 31), HaloID(3, 33), HaloID(3, 35), HaloID(3, 36), HaloID(3, 37), HaloID(3, 39), HaloID(3, 40),
#     HaloID(3, 41), HaloID(3, 43), HaloID(3, 45), HaloID(4, 1), HaloID(4, 2), HaloID(4, 4), HaloID(4, 6),
#     HaloID(4, 7), HaloID(4, 8), HaloID(4, 9), HaloID(4, 10), HaloID(4, 11), HaloID(4, 12), HaloID(4, 13),
#     HaloID(4, 16), HaloID(4, 17), HaloID(4, 19), HaloID(4, 21), HaloID(4, 24), HaloID(4, 25), HaloID(4, 27),
#     HaloID(4, 30), HaloID(4, 31), HaloID(4, 32), HaloID(4, 33), HaloID(4, 38), HaloID(4, 39), HaloID(4, 42),
#     HaloID(4, 44), HaloID(4, 47), HaloID(4, 48), HaloID(4, 49), HaloID(4, 50), HaloID(4, 52), HaloID(4, 56),
#     HaloID(4, 57), HaloID(4, 58), HaloID(5, 1), HaloID(5, 2), HaloID(5, 5), HaloID(5, 6), HaloID(5, 12),
#     HaloID(5, 22), HaloID(5, 27), HaloID(5, 30), HaloID(5, 48)
# ]

filename = data_path * "clusters_M200_gt_5e13.txt"
cluster_ids_Mgt5e4 = convert_halo_ids(filename)
store_data_path = data_path * "gamma_data/"

for cluster ∈ cluster_ids_Mgt5e4[3162:end]
    println(cluster)
    write_data(cluster, store_data_path)
end