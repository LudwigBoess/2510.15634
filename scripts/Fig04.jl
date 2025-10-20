include(joinpath(@__DIR__, "config.jl"))

using PyPlot, PyPlotUtility
using Printf
using SPHtoGrid
using GadgetIO, GadgetUnits
using Unitful, UnitfulAstro
using Cosmology


function plot_Fig04()

    clusters = ["Coma", "Virgo", "Perseus"]
    filenames = ["gamma_SB"]

    Ncols = length(clusters)
    Nrows = 2

    top_row = 1:Ncols
    bottom_row = (Ncols+1):2Ncols

    files = [map_path * "cic/$(cluster)_2rvir_000.$(filename).xy.fits"
             for cluster ∈ clusters, _ ∈ 1:2, filename ∈ filenames]

    im_cmap = ["plasma"]

    h = SnapshotHeader()
    h.h0 = 0.6777
    h.omega_0 = 0.307115
    h.omega_l = 0.692885
    h.z = 0.0
    GU = GadgetPhysical(h)

    # beam size: https://ui.adsabs.harvard.edu/abs/2009ApJ...697.1071A/abstract, Tab.1
    beam_size = 0.6u"°"

    z = [0.0231, 0.0038, 0.0179] #redshifts from SIMBAD

    smooth_sizes = Vector{Vector{Float64}}(undef, Ncols * Nrows)
    for i = 1:Ncols
        h.z = z[i]
        smooth_size = ustrip(arcmin_to_kpc(beam_size, h))
        smooth_sizes[i] =  [smooth_size, smooth_size]
        smooth_sizes[i+Ncols] = smooth_sizes[i]
    end
    annotate_smoothing = trues(Ncols*Nrows)
    annotate_smoothing[1:2:2Ncols] .= false
    smooth_file = trues(Ncols * Nrows)
    smooth_file[top_row] .= false

    cb_labels = ["Integrated " * L"\gamma" * "-ray Intensity  " * L"I_{\gamma, 0.2-300 \mathrm{ GeV}}" * " [erg s" * L"^{-1}" * " cm" * L"^{-2}" * "]"]

    annotate_time = trues(Nrows * Ncols)
    time_labels = ["" for _ = 1:Nrows*Ncols]
    time_labels[top_row] = clusters
    #time_labels[bottom_row] .= ""

    log_map = trues(Ncols * Nrows)

    vmin_arr = [1.e-13]
    vmax_arr = [1.e-8]

    transparent = false
    ticks_color = "k"

    plot_name = plot_path * "Fig04.pdf"

    plot_image_grid(Nrows, Ncols, files, im_cmap, cb_labels,
        vmin_arr, vmax_arr, plot_name, upscale=1.0,
        colorbar_mode="single",
        shift_colorbar_labels_inward=falses(Nrows),
        grid_direction="row",
        scale_label="1 Mpc",
        cb_label_offset=1.0,
        cutoffs=vmin_arr[1] .* ones(Ncols*Nrows),
        mask_bad=trues(Ncols * Nrows);
        smooth_sizes, annotate_smoothing, smooth_file,
        #par_arr, map_arr,
        time_labels, annotate_time,
        log_map,
        transparent, ticks_color
    )

end

plot_Fig04()
