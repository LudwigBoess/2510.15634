using GadgetIO, GadgetUnits
using SpectralCRsUtility
using Base.Threads
using ProgressMeter
using SPHtoGrid

#const global CR_factor = 1.60217733e-3 * 0.25

const global snap_base = "/e/ocean3/Local/3072/nonrad_mhd_crs_new/snapdir_000_z=0/snap_000"
const global sub_base = "/e/ocean3/Local/3072/nonrad_mhd_crs_new/groups_000_z=0/sub_000"


function get_gamma(data, h, Eγ;
    pmin=0.1, pmax=1.e5)

    center = [247.980, 245.480, 255.290] .* 1.e3

    GU = GadgetPhysical(h)

    # cr setup 
    Nbins = size(data["CRpN"], 1)
    par = CRMomentumDistributionConfig(pmin, pmax, Nbins)
    bounds = [par.pmin * 10.0^((i - 1) * par.bin_width) for i = 1:par.Nbins+1]

    Npart = size(data["CRpN"], 2)

    Sγ = Vector{Float64}(undef, Npart)

    @threads for i = 1:Npart

        norm = GU.CR_norm .* 10.0 .^ data["CRpN"][:, i]
        slope = Float64.(data["CRpS"][:, i])
        cut = Float64(data["CRpC"][i])

        pos = data["POS"][:, i] - center

        d = √(pos[1]^2 + pos[2]^2 + pos[3]^2) * GU.x_cgs
        #d = 19_000.0 * GU.x_cgs 

        nH = data["RHO"][i] * GU.rho_ncm3

        V = data["MASS"][i] / data["RHO"][i] * GU.x_cgs^3

        Sγ[i] = gamma_emissivity_pions(norm, slope, cut, bounds, nH, Eγ) * V / (4π * d^2)

    end

    return sum(Sγ)
end

function get_gamma_spectrum(data, h;
                    pmin=0.1, pmax=1.e5,
                    Nbins=20)

    
    Eγ = 10.0 .^ LinRange(-2, 4, Nbins)
    Lγ = Vector{Float64}(undef, Nbins)

    @showprogress for i = 1:length(Eγ)
        Lγ[i] = get_gamma(data, h, Eγ[i];
                            pmin, pmax)
    end

    return Eγ, Lγ
end

function get_gamma_flux(data, h, Eγ;
    pmin=0.1, pmax=1.e5)

    center = [247.980, 245.480, 255.290] .* 1.e3

    GU = GadgetPhysical(h)

    # cr setup 
    Nbins = size(data["CRpN"], 1)
    par = CRMomentumDistributionConfig(pmin, pmax, Nbins)
    bounds = [par.pmin * 10.0^((i - 1) * par.bin_width) for i = 1:par.Nbins+1]

    Npart = size(data["CRpN"], 2)

    Sγ = Vector{Float64}(undef, Npart)

    @threads for i = 1:Npart

        norm = GU.CR_norm .* 10.0 .^ data["CRpN"][:, i]
        slope = Float64.(data["CRpS"][:, i])
        cut = Float64(data["CRpC"][i])

        pos = data["POS"][:, i] - center

        d = √(pos[1]^2 + pos[2]^2 + pos[3]^2) * GU.x_cgs
        #d = 19_000.0 * 3.085678e21

        nH = data["RHO"][i] * GU.rho_ncm3

        V = data["MASS"][i] / data["RHO"][i] * GU.x_cgs^3

        Sγ[i] = gamma_flux_pions(norm, slope, cut, bounds, nH, V, d, Eγ_min=Eγ, Eγ_max=3.e3, N_integration_steps=20)

    end

    return sum(Sγ)
end

function get_gamma_flux_spectrum(data, h;
    pmin=0.1, pmax=1.e5,
    Nbins=20)


    Eγ = 10.0 .^ LinRange(log10(3e-2), log10(2e3), Nbins)
    Lγ = Vector{Float64}(undef, Nbins)

    @showprogress for i = 1:length(Eγ)
        Lγ[i] = get_gamma_flux(data, h, Eγ[i];
            pmin, pmax)
    end

    return Eγ, Lγ
end

function get_gamma_flux_PE(data, h, α_p, Eγ)


    center = [247.980, 245.480, 255.290] .* 1.e3

    GU = GadgetPhysical(h)

    Npart = length(data["RHO"])

    Sγ = Vector{Float64}(undef, Npart)

    @threads for i = 1:Npart

        pos = data["POS"][:, i] - center

        d = √(pos[1]^2 + pos[2]^2 + pos[3]^2) * GU.x_cgs
        #d = 19_000.0 * 3.085678e21

        rho_cgs = data["RHO"][i] * GU.rho_cgs
        T_K = data["U"][i] * GU.T_K 

        V = data["MASS"][i] / data["RHO"][i] * GU.x_cgs^3

        Sγ[i] = λγ_PE04(rho_cgs, T_K, α_p, Xcr=0.01, Eγ_π0_min=Eγ, Eγ_π0_max=2.e3) * V / (4π * d^2)

    end

    return sum(Sγ)
end

function get_gamma_flux_spectrum_PE(data, h, α_p;
    Nbins=20)


    Eγ = 10.0 .^ LinRange(log10(3e-2), log10(2e3), Nbins)
    Lγ = Vector{Float64}(undef, Nbins)

    @showprogress for i = 1:length(Eγ)
        Lγ[i] = get_gamma_flux_PE(data, h, α_p, Eγ[i])
    end

    return Eγ, Lγ
end


function get_gamma_flux_spectrum(halo_id::HaloID;
    pmin=0.1, pmax=1.e5,
    Nbins=20)

    blocks = ["POS", "MASS", "RHO", "U", "CRpP", "CRpN", "CRpS", "CRpC"]

    # only in r500
    gpos = read_subfind(sub_base * ".$(halo_id.file)", "GPOS")[:, halo_id.id]
    r200 = read_subfind(sub_base * ".$(halo_id.file)", "R200")[halo_id.id]

    sphere = GadgetSphere(gpos, r200)
    data = read_particles_in_geometry(snap_base, blocks, sphere, verbose=true)

    h = read_header(snap_base)

    return get_gamma_flux_spectrum(data, h;
        pmin, pmax,
        Nbins)
end

function get_gamma_flux_spectrum_PE04(halo_id::HaloID, α_p;
    Nbins=20)

    blocks = ["POS", "MASS", "RHO", "U"]

    # only in r500
    gpos = read_subfind(sub_base * ".$(halo_id.file)", "GPOS")[:, halo_id.id]
    r200 = read_subfind(sub_base * ".$(halo_id.file)", "R200")[halo_id.id]

    sphere = GadgetSphere(gpos, r200)
    data = read_particles_in_geometry(snap_base, blocks, sphere, verbose=true)

    h = read_header(snap_base)

    return get_gamma_flux_spectrum_PE(data, h, α_p;
        Nbins)
end

function get_gamma_spectrum(halo_id::HaloID;
                    pmin=0.1, pmax=1.e5,
                    Nbins=20)

    blocks = ["POS", "MASS", "RHO", "U", "CRpP", "CRpN", "CRpS", "CRpC"]

    # only in r500
    gpos = read_subfind(sub_base * ".$(halo_id.file)", "GPOS")[:, halo_id.id]
    r200 = read_subfind(sub_base * ".$(halo_id.file)", "R200")[halo_id.id]

    sphere = GadgetSphere(gpos, r200)
    data = read_particles_in_geometry(snap_base, blocks, sphere, verbose=true)

    h = read_header(snap_base)

    return get_gamma_spectrum(data, h;
        pmin, pmax,
        Nbins)
end

function write_gamma_spectrum(halo_id::HaloID, cluster_name::String)

    Eγ, Lγ = get_gamma_flux_spectrum(halo_id)

    filename = data_path * "gamma_data/gamma_$cluster_name.dat"
    f = open(filename, "w")
    write(f, length(Eγ))
    write(f, Eγ)
    write(f, Lγ)
    close(f)

end

function write_gamma_spectrum_PE04(halo_id::HaloID, cluster_name::String)

    for α_p in [2.2, 2.4, 2.5]
        Eγ, Lγ = get_gamma_flux_spectrum_PE04(halo_id, α_p)

        filename = data_path * "gamma_data/gamma_PE04_$cluster_name.$α_p.dat"
        f = open(filename, "w")
        write(f, length(Eγ))
        write(f, Eγ)
        write(f, Lγ)
        close(f)
    end
end

# Fornax
halo_names = ["Coma", "Virgo", "Perseus"]
#halo_ids = [HaloID(1, 11), HaloID(4, 44), HaloID(3, 9)] # old!
halo_ids = [HaloID(1, 9), HaloID(3, 39), HaloID(3, 12)]

for i = 1:3
    #write_gamma_spectrum(halo_ids[i], halo_names[i])
    write_gamma_spectrum_PE04(halo_ids[i], halo_names[i])
end