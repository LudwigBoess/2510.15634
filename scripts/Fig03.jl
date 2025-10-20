include(joinpath(@__DIR__, "config.jl"))

using GadgetIO, GadgetUnits
using SpectralCRsUtility
using Base.Threads
using Printf
using ProgressMeter
#using LoopVectorization
#using BenchmarkTools
using PyPlot, PyPlotUtility
using DelimitedFiles
using PyCall

const mp = 1.6726e-24
const cL = 2.9979e10
const erg2eV = 6.242e+11
conversion(p) = p * mp * cL^2 * erg2eV * 1.e-9

function plot_Fig03(r, CRpE, plot_name, transparent=false)

    axes_divider = pyimport("mpl_toolkits.axes_grid1.axes_divider")

    if transparent
        color = "w"
        cmap = PyPlot.cm.plasma_r
    else
        color="k"
        cmap = PyPlot.cm.plasma
    end

    Nbins = size(CRpE,2)

    par = CRMomentumDistributionConfig(0.1, 1.e5, Nbins)
    E_centers = conversion.(momentum_bin_centers(par))

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.LogNorm(vmin=1.e-2, vmax=2.0))
    sm.set_array([])


    fig = get_figure(1.0)
    plot_styling!(; color)
        ax = gca()
        axis_ticks_styling!(ax; color)
        ax.set_xscale("log")
        ax.set_yscale("log")
        xlabel("Proton Energy  " * L"E_p" * "  [GeV]")
        ylabel("Distr. Function  " * L"f\:(E_p)" * "  [arb. units]")
        ax.set_xlim([1.e-1, 1.e4])
        ax.set_ylim([5.e-20, 3.e-11])

        axvline(1.22, linestyle="--", color=color, alpha=0.7)

        for i = Nbins:-1:1 
            plot(E_centers, CRpE[i,:] ./ E_centers, color=sm.to_rgba(r[i]))
        end

        α = 2.0
        f0 = 1.e-13
        E1 = 10.0
        E2 = 2.e2
        plot([E1, E2], [f0, f0 * (E1 / E2)^(α)], color=color, linestyle="--")

        text(30.0, 5.e-15, L"\alpha = -2", rotation=-49)
        text(1.5, 5.e-19, L"T_{p,\mathrm{thr.}}")

    ax_divider = axes_divider.make_axes_locatable(ax)
    cax = ax_divider.append_axes("right", size = "7%", pad = "2%")
    cb = colorbar(sm, cax=cax, location="right")#, anchor=(0.1, 0.5))
    cb_ticks_styling!(cb; color)
    cb.set_label("Radius  " * L"r/r_{200}")

    savefig(plot_name, bbox_inches="tight", dpi=400; transparent)
    close(fig)
end

data = readdlm(data_path * "CRpE_coma.txt")
r_bins = data[:,1]
CRpE_bins = data[:,2:end]
r_plot = 10.0.^r_bins
CRpE_plot = 10.0.^CRpE_bins
plot_name = plot_path * "Fig03.pdf"
println(r_plot)

plot_Fig03(r_plot, CRpE_plot, plot_name, false)

