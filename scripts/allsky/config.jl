# path to the package envirnoment
const global dir_path = joinpath(@__DIR__, "..")

# activate the environment
using Pkg;
Pkg.activate(dir_path);

# simulation settings
const global snap_base = "/path/to/snapshot"

# general map settings
const global Nside = 1024
const global Nfiles = 2048
const global map_path = "/path/to/maps/allsky/"
const global viewpoint = "slow_1"
const global kernel = WendlandC4(2)

const global h = GadgetIO.read_header(snap_base)
const global GU = GadgetPhysical(h)
const global kpc2cm = 3.085678e21


# SLOW 1 paper
const global center_comov = [247.980, 245.480, 255.290] .* 1.e3
const global center = center_comov .* GU.x_physical
global const radius_limits = [5_000.0, 240_000.0 * GU.x_physical]
#const global radius_limits = [40_000.0, 300_000.0] # for CReE

"""
    find_outside_shell(pos)

Find particles outside shell with given radius limits.
"""
function find_outside_shell(pos)
    # calculate radii of all particles
    Δx = Vector{Float64}(undef, size(pos, 2))

    @threads for i = 1:length(Δx)
        Δx[i] = 0.0
        for dim = 1:3
            Δx[i] += (pos[dim, i] - center[dim])^2
        end
        Δx[i] = sqrt(Δx[i])
    end

    @. ((Δx < radius_limits[1]) || (Δx > radius_limits[2]))
end

"""
    set_rest_to_zero(pos, quantitiy)

Set all particles outside shell to zero.
"""
function set_rest_to_zero(pos, quantitiy)

    sel = find_outside_shell(pos)

    quantitiy[sel] .= 0.0

    return quantitiy
end