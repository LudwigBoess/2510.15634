function CRpP_maps_of_subfile(subfile)

    println("CRpP: subfile $subfile running on $(nthreads()) threads")
    flush(stdout); flush(stderr)

    hsml = read_block(snap_base * ".$subfile", "HSML", parttype=0) .* GU.x_physical
    rho  = read_block(snap_base * ".$subfile", "RHO", parttype=0)  .* GU.rho_physical
    m   = read_block(snap_base * ".$subfile", "MASS", parttype=0) .* GU.m_physical  

    CRpP = read_block(snap_base * ".$subfile", "CRpP", parttype=0) .* GU.P_CR_cgs

    pos  = read_block(snap_base * ".$subfile", "POS", parttype=0) .* GU.x_physical
    map = healpix_map(pos, hsml, m, rho, CRpP, rho, show_progress=false, calc_mean=false;
        center, kernel, Nside)

    pos  = hsml = rho = m = nothing
    CRpP = nothing
    GC.gc()

    return map
end


function Xcr_maps_of_subfile(subfile)

    println("Xcr: subfile $subfile running on $(nthreads()) threads")
    flush(stdout); flush(stderr)

    hsml    = read_block(snap_base * ".$subfile", "HSML", parttype=0) .* GU.x_physical
    rho = read_block(snap_base * ".$subfile", "RHO", parttype=0)

    CRpP = read_block(snap_base * ".$subfile", "CRpP", parttype=0)

    U = read_block(snap_base * ".$subfile", "U", parttype=0)

    γ_m1 = 5/3 - 1
    Xcr = Vector{Float64}(undef, length(hsml))

    @threads for i = 1:length(hsml)
        Xcr[i]  = CRpP[i] / (γ_m1 * rho[i] * U[i])
    end

    pos  = read_block(snap_base * ".$subfile", "POS", parttype=0) .* GU.x_physical
    map = healpix_map(pos, hsml, m, rho, Xcr, rho, show_progress=false, calc_mean=false;
                        center, kernel, Nside)

    pos  = hsml = rho = m = nothing
    CRpP = nothing
    Xcr  = nothing
    U    = nothing
    GC.gc()
    
    return map
end

function Xcr_mean_maps_of_subfile(subfile)

    println("Xcr: subfile $subfile running on $(nthreads()) threads")
    flush(stdout); flush(stderr)

    hsml    = read_block(snap_base * ".$subfile", "HSML", parttype=0) .* GU.x_physical
    rho = read_block(snap_base * ".$subfile", "RHO", parttype=0)

    CRpP = read_block(snap_base * ".$subfile", "CRpP", parttype=0)

    U = read_block(snap_base * ".$subfile", "U", parttype=0)

    γ_m1 = 5/3 - 1
    Xcr = Vector{Float64}(undef, length(hsml))

    @threads for i = 1:length(hsml)
        Xcr[i]  = CRpP[i] / (γ_m1 * rho[i] * U[i])
    end

    pos  = read_block(snap_base * ".$subfile", "POS", parttype=0) .* GU.x_physical
    map = healpix_map(pos, hsml, m, rho, Xcr, rho, show_progress=false, calc_mean=true;
                        center, kernel, Nside)

    pos  = hsml = rho = m = nothing
    CRpP = nothing
    Xcr  = nothing
    U    = nothing
    GC.gc()
    
    return map
end