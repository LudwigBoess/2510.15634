using Healpix
using SkyCoords
using ProgressMeter
#using PyPlot, PyPlotUtility
using Printf
using Base.Threads


function supergalactic_angle(res, i)

    b, l = pix2angRing(res, i)
    gal = GalCoords(l, colat2lat(b))
    supergal = convert(SuperGalCoords, gal)

    return (π / 2 - supergal.b), supergal.l
end


function rotate_healpix_map(filename, out_file)

    m = Healpix.readMapFromFITS(filename, 1, Float64)
    Nside = npix2nside(length(m[:]))

    res = Healpix.Resolution(Nside)

    m_supergal = HealpixMap{Float64,RingOrder}(Nside)

    P = Progress(nside2npix(Nside))

    for i ∈ 1:nside2npix(Nside)
        θ, ϕ = supergalactic_angle(res, i)
        m_supergal[i] = interpolate(m, θ, ϕ)

        next!(P)
    end

    if isfile(out_file)
        rm(out_file)
    end
    saveToFITS(m_supergal, out_file)
    m_supergal = nothing
end


function rotate_map(map_path, filename, out_file)
    rotate_healpix_map(map_path * filename, map_path * out_file)
end

function rotate_all_files(map_path)

    files = readdir(map_path)


    for filename ∈ files

        println(filename)

        split_filename = split(filename, ".")

        out_file = split_filename[1] * "_gal.fits"
        if !isfile(out_file)
            rotate_map(map_path, filename, out_file)
        end
    end
end

map_path = "/path/to/maps/allsky/"
filenames = ["allsky_gammaSB_slow_1.fits", "allsky_CRpP_slow_1.fits"]

for fi in filenames 

    println(fi)
    split_filename = split(fi, ".")

    out_file = split_filename[1] * "_gal.fits"
    #if !isfile(out_file)
    rotate_map(map_path, fi, out_file)
    println("done!")
end
