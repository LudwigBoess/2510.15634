include(joinpath(@__DIR__, "config.jl"))

using PyPlot, PyPlotUtility
using Printf
using GadgetIO, GadgetUnits
using DiffusiveShockAccelerationModels

function read_Fkin(filename)
    f = open(filename, "r")
    Nbins = read(f, Int64)
    bin_lim = read!(f, Vector{Float64}(undef,  2))
    mach_count = read!(f, Vector{Int64}(undef, Nbins))
    Fkin = read!(f, Vector{Float64}(undef, Nbins))
    Fkin_weighted = read!(f, Vector{Float64}(undef, Nbins))
    close(f)

    bin_width = (bin_lim[2] - bin_lim[1]) / Nbins
    mach_centers = [10.0^bin_lim[1] * 10.0^((i - 0.5) * bin_width) for i = 1:Nbins]

    dM = [10.0^(i * bin_width) - 10.0^((i - 1) * bin_width) for i = 1:Nbins]

    return mach_centers, dM, Fkin, Fkin_weighted
end


function plot_Fig06(mach, dM, Fkin, Fkin_weighted, plot_name)

    lw = 3

    fig = get_figure()
    plot_styling!()

    ax = gca()
    axis_ticks_styling!(ax)
    ax.set_xlim([1.0, 1.e3])
    ax.set_ylim([1.e-6, 1.e3])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Sonic Mach Number " * L"\mathcal{M}_s")
    #ax.set_ylabel("Energy Flux " * L"F" * "  [erg s" * L"^{-1}" * " cm" * L"^{-2}" * "]")
    ax.set_ylabel("Energy Flux " * L" \frac{\mathrm{d}F}{\mathrm{d}\mathcal{M}_s}" * "  [erg s" * L"^{-1}" * " cm" * L"^{-2}" * "]")

    
    plot(mach, Fkin, color="purple", lw=lw, label=L"F_\mathrm{kin}")
    plot(mach, Fkin_weighted, color="purple", linestyle=ls("__"), lw=lw, label=L"\eta(\theta_B) \times F_\mathrm{kin}")
    η = @. η_Ms_reacc(Ryu19(), mach)
    plot(mach, η .* Fkin_weighted, color="purple", linestyle="--", lw=lw, label=L"\eta(\theta_B) \times \eta(\mathcal{M}_s, X_\mathrm{cr}) \times F_\mathrm{kin}")
    η = @. η_Ms_acc(Ryu19(), mach)
    plot(mach, η .* Fkin_weighted, color="purple", linestyle=":", lw=lw, label=L"\eta(\theta_B) \times \eta(\mathcal{M}_s, 0) \times F_\mathrm{kin}")

    legend(frameon=false, loc="upper right")
    legend(frameon=false, bbox_to_anchor=(0.5, 1.0), ncol=2, loc="lower center")

    savefig(plot_name, bbox_inches="tight")
    close(fig)
end

filename = data_path * "Fkin.dat"
mach, dM, Fkin, Fkin_weighted = read_Fkin(filename)
plot_name = plot_path * "Fig06.pdf"
plot_Fig06(mach, dM, Fkin, Fkin_weighted, plot_name)
