function get_gamma(data, calc_flux=false)

    Npart = length(data["CRpC"])

    d_min = radius_limits[1] * GU.x_cgs
    d_max = radius_limits[2] * GU.x_cgs

    Fγ = Vector{Float64}(undef, Npart)

    # cr setup 
    Nbins = size(data["CRpN"], 1)
    par = CRMomentumDistributionConfig(0.1, 1.e5, Nbins)
    # construct boundaries 
    bounds = momentum_bin_boundaries(par)

    @threads for i ∈ eachindex(Fγ)

        d = 0.0
        @inbounds for dim = 1:3
            d += (data["POS"][dim, i] - center_comov[dim])^2
        end

        d = sqrt(d) * GU.x_cgs

        if d_min <= d <= d_max
            
            norm = GU.CR_norm .* 10.0 .^ data["CRpN"][:, i]

            slope = Float64.(data["CRpS"][:, i])
            cut = Float64(data["CRpC"][i])

            nH = data["RHO"][i] * GU.rho_ncm3

            V = 1.0

            if calc_flux
                Fγ[i] = gamma_flux_pions(norm, slope, cut, bounds, nH, V, d, N_integration_steps=20,
                    Eγ_min=0.5, Eγ_max=200.0)
            else
                Fγ[i] = gamma_luminosity_pions(norm, slope, cut, bounds, nH, V, N_integration_steps=20,
                    Eγ_min=0.5, Eγ_max=200.0)
            end

        else
            Fγ[i] = 0.0
        end

    end

    Fγ
end

function gamma_maps_of_subfile(subfile)

    println("gamma: subfile $subfile running on $(nthreads()) threads")
    flush(stdout)
    flush(stderr)

    hsml = read_block(snap_base * ".$subfile", "HSML", parttype=0) .* GU.x_physical
    rho = read_block(snap_base * ".$subfile", "RHO", parttype=0) .* GU.rho_physical
    m = read_block(snap_base * ".$subfile", "MASS", parttype=0) .* GU.m_physical

    weights = part_weight_physical(length(hsml), GU.x_cgs)

    data = Dict(block => read_block(snap_base * ".$subfile", block, parttype=0)
                for block ∈ ["POS", "MASS", "RHO", "CRpN", "CRpS", "CRpC"])

    h = read_header(snap_base * ".$subfile")
    Fγ = get_gamma(data)

    println("\tgamma done!\n\tmaximum = $(maximum(Fγ)) erg/s/cm^3\n\tsum = $(sum(Fγ)) erg/s/cm^3")
    flush(stdout)
    flush(stderr)

    pos = read_block(snap_base * ".$subfile", "POS", parttype=0) .* GU.x_physical

    map = healpix_map(pos, hsml, m, rho, Fγ, weights, show_progress=false;
        center, kernel, Nside, radius_limits)

    pos = hsml = rho = m = weights = nothing
    Fγ = nothing
    data = nothing
    GC.gc()

    return map

end

