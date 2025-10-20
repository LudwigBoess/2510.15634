"""
    Reproduce Fig. 05 from Ha+2020: 
    https://ui.adsabs.harvard.edu/abs/2020ApJ...892...86H/abstract

    Based on data by Ackermann+14:
    https://iopscience.iop.org/article/10.1088/0004-637X/787/1/18
    
    Gamma energy range: 500MeV - 200GeV
"""


include(joinpath(@__DIR__, "config.jl"))

@info "loading packages"
using GadgetIO, GadgetUnits
using SpectralCRsUtility
using Statistics
using ProgressMeter
using Base.Threads
using Statistics
using DelimitedFiles
@info "loading PyPlot"
using PyPlot, PyPlotUtility
@info "PyPlot done"
using Unitful, UnitfulAstro
using Cosmology
@info "all done"


function convert_ackermann(data_path)

    ackermann = readdlm(data_path * "gamma_data/ackermann.txt", ',', comments=true)
    m = ackermann[:,2] .* 0.7 ./ 0.677  # compensate for different little h

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

function plot_Fig02(data_path, plot_name, transparent=false)

    if transparent
        color = "w"
        dot_color = "#8f0ec7"
        x_color = "#c70eb5"
    else
        color = "k"
        dot_color = "#8f0ec7"
        x_color = "#c70eb5"
    end


    # read
    ackermann = readdlm(data_path * "gamma_data/ackermann.txt", ',', comments=true)
    m, Lγ = convert_ackermann(data_path)


    # read simulation data
    sim = readdlm(data_path * "gamma_data/gamma_R200.txt")

    # bin in mass
    Nbins = 20
    mass_bins = LinRange(log10(5.e13), log10(5.e15), Nbins)
    mass_plot_bins = 10.0.^mass_bins

    # bin mean in log space
    mass_log10 = log10.(sim[:, 2])
    Xcr_mean = 10.0.^bin_1D(mass_log10, [mass_bins[1], mass_bins[end]], log10.(sim[:,3]), calc_mean=true, calc_sigma=false, Nbins=Nbins)[2]
    Xcr_gt0_mean = 10.0.^bin_1D(mass_log10, [mass_bins[1], mass_bins[end]], log10.(sim[:,4]), calc_mean=true, calc_sigma = false, Nbins=Nbins)[2]
    Fgamma_mean = 10.0.^bin_1D(mass_log10, [mass_bins[1], mass_bins[end]], log10.(sim[:,5]), calc_mean=true, calc_sigma = false,  Nbins=Nbins)[2]
    Fgamma_gt0_mean = 10.0.^bin_1D(mass_log10, [mass_bins[1], mass_bins[end]], log10.(sim[:,5].* sim[:, 4] ./ sim[:, 3]), calc_mean=true, calc_sigma = false, Nbins=Nbins)[2]
    Lgamma_mean = 10.0.^bin_1D(mass_log10, [mass_bins[1], mass_bins[end]], log10.(sim[:,6]), calc_mean=true, calc_sigma = false, Nbins=Nbins)[2]
    Lgamma_gt0_mean = 10.0.^bin_1D(mass_log10, [mass_bins[1], mass_bins[end]], log10.(sim[:,6].* sim[:, 4] ./ sim[:, 3]), calc_mean=true, calc_sigma = false, Nbins=Nbins)[2]
    
    # set zeros to NaN for plotting
    for dummy in [Xcr_mean, Xcr_gt0_mean, Fgamma_mean, Fgamma_gt0_mean, Lgamma_mean, Lgamma_gt0_mean]
        dummy[dummy .== 1.0] .= NaN
    end

    x_pixels = 1400
    axis_label_font_size = 8
    legend_font_size = 6
    ideal_alpha = 0.3
    fig = get_figure(0.3; x_pixels)
    plot_styling!(x_pixels; axis_label_font_size, legend_font_size)

    gs = plt.GridSpec(3, 1, figure=fig, hspace=0.01)

    subplot(get_gs(gs, 0, 0))

    ax = gca()
    axis_ticks_styling!(ax)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([1.e13, 5.e15])
    ax.set_ylim([5.e-7, 1.e-1])
    ax.set_xticklabels([])
    ax.set_ylabel("Pressure ratio\n" * L"X_{cr}")
    ax.yaxis.set_label_coords(-0.2, 0.5)

    errorbar(m, ackermann[:, 3],
        xerr=0.1 .* m, yerr=0.3 .* ackermann[:, 3],
        fmt=",", uplims=trues(length(ackermann[:, 1])),
        color=color, alpha=1.0)
    scatter(sim[:, 2], sim[:, 3], marker="o", rasterized=true, alpha=ideal_alpha, s=15, linewidths=0.0, color=dot_color)
    scatter(sim[:, 2], sim[:, 4], marker="x", rasterized=true, alpha=ideal_alpha, s=9, color=x_color)
    scatter([0.0], [0.0], marker="o", color=dot_color, label=L"\mathcal{Q}")
    scatter([0.0], [0.0], marker="x", color=x_color, label=L"\chi \times \mathcal{Q}")
    plot(mass_plot_bins, Xcr_mean, color="white", alpha=0.5, lw=5)
    plot(mass_plot_bins, Xcr_mean, color=dot_color, label=L"<\log_{10}(\mathcal{Q})>", lw=3)
    plot(mass_plot_bins, Xcr_gt0_mean, color="white", alpha=0.5, lw=5)
    plot(mass_plot_bins, Xcr_gt0_mean, linestyle="dashed", color=x_color, label=L"<\log_{10}(\chi \times \mathcal{Q})", lw=3)

    scatter([0.0], [0.0], c="k", marker=L"\downarrow", s=80, label="Ackermann+14")
    legend(frameon=false, loc="lower center", bbox_to_anchor=(0.45, 1.05), ncols=3)


    subplot(get_gs(gs, 1, 0))

    ax = gca()
    axis_ticks_styling!(ax)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([1.e13, 5.e15])
    ax.set_xticklabels([])
    ax.set_ylim([2.e-17, 3.e-9])
    ax.set_ylabel("γ-ray Flux\n" * L"F_{\gamma \in [0.5-200] \: \mathrm{GeV}}" * " [γ cm" * L"^{-2}" * "s" * L"^{-1}" * "]")
    ax.yaxis.set_label_coords(-0.2, 0.5)

    errorbar(m, ackermann[:, 5],
        xerr=0.1 .* m, yerr=0.3 .* ackermann[:, 5],
        fmt=",", uplims=trues(length(ackermann[:, 1])),
        color=color, alpha=1.0)
    scatter(sim[:, 2], sim[:, 5], marker="o", alpha=ideal_alpha, s=15, linewidths=0.0, color=dot_color, zorder=2, rasterized=true)
    scatter(sim[:, 2], sim[:, 5] .* sim[:, 4] ./ sim[:, 3], marker="x", alpha=ideal_alpha, s=9, color=x_color, zorder=2, rasterized=true)
    plot(mass_plot_bins, Fgamma_mean, color="white", alpha=0.5, lw=5)
    plot(mass_plot_bins, Fgamma_mean, color=dot_color, label=L"<Q>", lw=3)
    plot(mass_plot_bins, Fgamma_gt0_mean, color="white", alpha=0.5, lw=5)
    plot(mass_plot_bins, Fgamma_gt0_mean, linestyle="dashed", color=x_color, label=L"<Q>", lw=3)

    subplot(get_gs(gs, 2, 0))

    ax = gca()
    axis_ticks_styling!(ax; color)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([1.e13, 5.e15])
    ax.set_xlabel("Mass  " * L"M_{200}" * " [" * L"M_\odot" * "]")
    ax.set_ylim([2.e38, 2.e46])
    ax.set_ylabel("γ-ray Luminosity\n" * L"L_{\gamma \in [0.5-200] \: \mathrm{GeV}}" * " [γ s" * L"^{-1}" * "]")
    ax.yaxis.set_label_coords(-0.2, 0.5)

    errorbar(m, Lγ,
        xerr=0.1 .* m, yerr=0.3 .* Lγ,
        fmt=",", uplims=trues(length(Lγ)),
        color=color, alpha=1.0)
    scatter(sim[:, 2], sim[:, 6], marker="o", alpha=ideal_alpha, s=15, linewidths=0.0, color=dot_color, zorder=2, rasterized=true)
    scatter(sim[:, 2], sim[:, 6] .* sim[:, 4] ./ sim[:, 3], marker="x", alpha=ideal_alpha, s=9, color=x_color, zorder=2, rasterized=true)
    plot(mass_plot_bins, Lgamma_mean, color="white", alpha=0.5, lw=5)
    plot(mass_plot_bins, Lgamma_mean, color=dot_color, label=L"<Q>", lw=3)
    plot(mass_plot_bins, Lgamma_gt0_mean, color="white", alpha=0.5, lw=5)
    plot(mass_plot_bins, Lgamma_gt0_mean, linestyle="dashed", color=x_color, label=L"<Q>", lw=3)

    @info "saving"
    savefig(plot_name, bbox_inches="tight", dpi=400)
    close(fig)
end

plot_name = plot_path * "Fig02.pdf"
plot_Fig02(data_path, plot_name, false)

