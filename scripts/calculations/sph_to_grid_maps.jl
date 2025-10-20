using Distributed
addprocs(4)

using GadgetIO, GadgetUnits
using SPHKernels, SPHtoGrid
using Printf
using ProgressMeter
using SpectralCRsUtility
using Base
using LinearAlgebra
using Base.Threads
using Rotations

global const snap_base = "/e/ocean3/Local/3072/nonrad_mhd_crs_new/snapdir_000_z=0/snap_000"
global const slow_center = [247.980, 245.480, 255.290] .* 1.e3


function get_gamma(data, h)

    GU = GadgetPhysical(h)

    Npart = length(data["CRpC"])
    jγ = Vector{Float64}(undef, Npart)

    # cr setup 
    Nbins = size(data["CRpN"], 1)
    par = CRMomentumDistributionConfig(0.1, 1.e5, Nbins)
    # construct boundaries 
    bounds = momentum_bin_boundaries(par)

    p = Progress(Npart)

    @threads for i = 1:Npart

        norm = GU.CR_norm .* 10.0 .^ data["CRpN"][:, i]

        slope = Float64.(data["CRpS"][:, i])
        cut = Float64(data["CRpC"][i])

        nH = data["RHO"][i] * GU.rho_ncm3

        # [erg / s / cm^3]
        jγ[i] = gamma_luminosity_pions(norm, slope, cut, bounds, nH, 1.0)

        next!(p)
    end

    jγ
end

function make_maps(snap, cluster, gpos, rvir)

    blocks = ["POS",  "HSML", "RHO", "U", "MASS",
                "CRpN", "CRpS", "CRpC"
            ]

    @info "reading data"

    # default
    image_path = map_path * "cic/$(cluster)_2rvir_$(@sprintf("%03i", snap))."
    h = read_header(snap_base)
    cube = GadgetCube(gpos, rvir)
    data = read_particles_in_geometry(snap_base, blocks, cube, use_keys=true)

    @info "done"
    # select kernel
    kernel = WendlandC4(Float64, 2)

    GU = GadgetPhysical(h)

    # project along LOS
    a = [1, 1, 1] ./ sqrt(3)
    b = gpos - slow_center
    b /= norm(b)

    I = one(RotMatrix{3, Float64})
    v = a × b
    c = a ⋅ b
    V = [    0  -v[3]  v[2] 
           v[3]    0  -v[1]
          -v[2]  v[1]    0 ]
    R = I + V + V^2 * (1 / (1 + c))

    # convert to physical code units for mapping
    pos = R * (data["POS"] .- gpos) .* GU.x_physical
    hsml = data["HSML"] .* GU.x_physical
    rho = data["RHO"] .* GU.rho_physical
    mass = data["MASS"] .* GU.m_physical

    m_cgs = data["MASS"] .* GU.m_cgs

    xy_size = 2rvir
    z_size = 2rvir

    Npart = length(m_cgs)

    # define mapping parameters
    param = mappingParameters(center=zeros(3),
        x_size=xy_size * GU.x_physical,
        y_size=xy_size * GU.x_physical,
        z_size=z_size * GU.x_physical,
        Npixels=1024)

    # Synch
    @info "gamma"
    image_prefix = image_path * "gamma_SB"

    jγ = get_gamma(data, h)
    weights = part_weight_physical(length(jγ), param)
    units = "erg/s/cm^2"
    reduce_image = false

    map_it(pos, hsml, mass, rho, jγ, weights,
        parallel=true,
        projection="xy";
        reduce_image, units,
        kernel, snap, param, image_prefix)


end

# coma
gpos = [245040.45, 327781.84, 246168.69]
rvir = 2177.5625
cluster = "Coma"
make_maps(0, cluster, gpos, 1.5rvir)


# virgo
gpos = [244450.58, 255851.78, 253668.19]
rvir = 1711.3118
cluster = "Virgo"
make_maps(0, cluster, gpos, 1.5rvir)


# perseus
gpos = [307184.78, 247627.08, 230736.48]
rvir = 1813.042
cluster = "Perseus"
make_maps(0, cluster, gpos, 1.5rvir)
