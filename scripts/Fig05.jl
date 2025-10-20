include(joinpath(@__DIR__, "config.jl"))

using GadgetIO, GadgetUnits
using SpectralCRsUtility
using PyPlot, PyPlotUtility
using Base.Threads
using ProgressMeter
using ColorSchemes
using SPHtoGrid


function read_gamma_spectrum(cluster_name)
    filename = data_path * "gamma_data/gamma_$cluster_name.dat"
    f = open(filename, "r")
    Nbins = read(f, Int64)
    Eγ = read!(f, Vector{Float64}(undef, Nbins))
    Lγ = read!(f, Vector{Float64}(undef, Nbins))
    close(f)

    return Eγ, Lγ
end

function read_gamma_spectrum_PE(cluster_name, slope)
    filename = data_path * "gamma_data/gamma_PE04_$cluster_name.$slope.dat"
    f = open(filename, "r")
    Nbins = read(f, Int64)
    Eγ = read!(f, Vector{Float64}(undef, Nbins))
    Lγ = read!(f, Vector{Float64}(undef, Nbins))
    close(f)

    return Eγ, Lγ
end

function get_spectral_slope(Eγ, Lγ)

    Nbins = length(Eγ)
    α = Vector{Float64}(undef, Nbins-1)
    Eγ_center = Vector{Float64}(undef, Nbins - 1)

    for i = 1:Nbins-1
        Eγ_center[i] = 10.0^(0.5 * (log10(Eγ[i]) + log10(Eγ[i+1])))
        α[i] = abs((log10(Lγ[i+1]) - log10(Lγ[i])) / (log10(Eγ[i+1]) - log10(Eγ[i])))
    end

    return Eγ_center, α
end

const perseus_limits = [ 630.0 3.22e-13
                        1000.0 1.38e-13
                        1600.0 1.18e-13
                        2500.0 0.87e-13 ]

const perseus_aleksic = [100.0 6.55e-12
                        130.0 6.21e-12
                        160.0 6.17e-12
                        200.0 5.49e-12
                        250.0 4.59e-12
                        320.0 3.36e-12
                        400.0 1.83e-12
                        500.0 1.39e-12
                        630.0 0.72e-12
                        800.0 0.65e-12
                        1000.0 0.472e-12]



function plot_Fig05(names, plot_name)

    #colors = ["purple", "teal", "darkblue", "blue"]
    #colors = [1, 3, 5, 7]
    colors = collect(1:4)
    colors = collect(9:-1:6)

    x_pixels = 1000
    fig = get_figure(3.0; x_pixels)
    plot_styling!(x_pixels)

    gs = plt.GridSpec(2, 4, figure=fig, height_ratios=[1, 0.4], hspace=0.05, wspace=0.05)  

    for i = 1:length(names)

        subplot(get_gs(gs, 0, i-1))
        ax1 = gca()
        axis_ticks_styling!(ax1, tick_label_size=22)
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.set_ylim([4.e-15, 5.e-7])
        ax1.set_xticklabels([])
        ax1.set_xlim([3.e-2, 2.e3])
        locmaj = matplotlib.ticker.LogLocator(base=10, numticks=12)
        ax1.yaxis.set_major_locator(locmaj)
        if i == 1
            ax1.set_ylabel(L"F_\gamma(>E_\gamma)" * " [" * L"\gamma" * " s" * L"^{-1}" * " cm" * L"^{-2}" * "]")
            ax1.yaxis.set_label_coords(-0.14, 0.5)
        else
            ax1.set_yticklabels([])
        end


        subplot(get_gs(gs, 1, i - 1))
        ax2 = gca()
        axis_ticks_styling!(ax2, tick_label_size=22)
        ax2.set_xscale("log")
        #ax.set_yscale("log")
        ax2.set_ylim([0.0, 3.0])
        ax2.set_xlim([3.e-2, 2.e3])
        ax2.set_xlabel("Photon Energy  " * L"E_\gamma" * " [GeV]")
        locmaj = matplotlib.ticker.LogLocator(base=10, numticks=12)
        ax2.xaxis.set_major_locator(locmaj)
        if i == 1
            ax2.set_ylabel(#"Spectral Slope  " * 
                            L"\alpha_\gamma")
            ax2.yaxis.set_label_coords(-0.14, 0.5)
        else
            ax2.set_yticklabels([])
        end

        # ax1.axvspan(0.1, 300.0, color="k", alpha=0.2)
        # ax1.text(1.0, 1.e-8, "Fermi-LAT")

        Eγ, Lγ = read_gamma_spectrum(names[i])
        ax1.plot(Eγ, Lγ, color="purple", lw=5)
        Eγ_center, α = get_spectral_slope(Eγ, Lγ)
        ax2.plot(Eγ_center, α, color="purple", lw=5)

        Eγ, Lγ = read_gamma_spectrum_PE(names[i], 2.5)
        ax1.plot(Eγ[1:end-1], Lγ[1:end-1], color="purple", linestyle="--", lw=2)
        Eγ_center, α = get_spectral_slope(Eγ[1:end-1], Lγ[1:end-1])
        ax2.plot(Eγ_center, α, color="purple", linestyle="--", lw=2)

        Eγ, Lγ = read_gamma_spectrum_PE(names[i], 2.4)
        ax1.plot(Eγ[1:end-1], Lγ[1:end-1], color="purple", linestyle=":", lw=2)
        Eγ_center, α = get_spectral_slope(Eγ[1:end-1], Lγ[1:end-1])
        ax2.plot(Eγ_center, α, color="purple", linestyle=":", lw=2)

        Eγ, Lγ = read_gamma_spectrum_PE(names[i], 2.2)
        ax1.plot(Eγ[1:end-1], Lγ[1:end-1], color="purple", linestyle=ls(".-"), lw=2)
        Eγ_center, α = get_spectral_slope(Eγ[1:end-1], Lγ[1:end-1])
        ax2.plot(Eγ_center, α, color="purple", linestyle=ls(".-"), lw=2)

        #ax1.text(150.0, 5.e-9, names[i])
        ax1.set_title(names[i])

        if names[i] == "Coma"

            ax1.errorbar([0.2], [3.08e-9],
                xerr=0.1 .* [0.2], yerr=0.3 .* [3.08e-9],
                color="k",
                fmt=",", uplims=trues(1))

            ax1.scatter([0.0], [0.0], c="k", marker=L"\downarrow", s=80, label="Xi+18")

        end

        if names[i] == "Virgo"
            
            ax1.errorbar([0.1], [15.0e-9],
                xerr=0.1 .* [0.1], yerr=0.3 .* [15.0e-9],
                color="k",
                fmt=",", uplims=trues(1))

            ax1.errorbar([1.e3], [1.50e-12],
                xerr=0.1 .* [1.e3], yerr=0.3 .* [1.50e-12],
                color="gray",
                fmt=",", uplims=trues(1))

            ax1.scatter([0.0], [0.0], c="k", marker=L"\downarrow", s=80, label="Ackermann+15")
            ax1.scatter([0.0], [0.0], c="gray", marker=L"\downarrow", s=80, label="H.E.S.S.+23")
        end


        if names[i] == "Perseus"
            ax1.errorbar(perseus_aleksic[:, 1], perseus_aleksic[:, 2],
                xerr=0.1 .* perseus_aleksic[:, 1], 
                yerr=0.3 .* perseus_aleksic[:, 2],
                        color="k", 
                        fmt=",", uplims=trues(length(perseus_aleksic[:, 2])))

            ax1.errorbar(perseus_limits[:, 1], perseus_limits[:, 2],
                xerr=0.1 .* perseus_limits[:, 1], yerr=0.3 .* perseus_limits[:, 2],
                color="gray",
                fmt=",", uplims=trues(length(perseus_limits[:, 2])))
            
            ax1.scatter([0.0], [0.0], c="k", marker=L"\downarrow", s=80, label="Aleksic+11")
            ax1.scatter([0.0], [0.0], c="gray", marker=L"\downarrow", s=80, label="Aleksic+12")


        end

        if names[i] == "Fornax"
            ax1.plot(pinzke_data[:, 1], pinzke_data[:, 2], color="gray", linestyle="--", lw=4, label="Pinzke+11")
            Eγ_center, α = get_spectral_slope(pinzke_data[:, 1], pinzke_data[:, 2])
            ax2.plot(Eγ_center, α, color="gray", linestyle="--", lw=4)

            ax1.errorbar([0.1], [0.54e-9], xerr=0.1 .* [0.1], yerr=0.3 .* [0.54e-9],
                        color="k", fmt=",", uplims=[true])
            ax1.scatter([0.0], [0.0], c="k", marker=L"\downarrow", s=80, label="Ando+12")
        end

        ax1.legend(frameon=false)#, bbox_to_anchor=(0.5, 1.0))

    end

    subplot(get_gs(gs, 1, 1))
    ax1 = gca()
    ax1.plot([0.0], [0.0], color="purple", lw=5,label="This work")
    ax1.plot([0.0], [0.0], color="purple", lw=2, linestyle="--", label="Pfrommer&Ensslin (2004) " * L"\alpha_p = 2.5")
    ax1.plot([0.0], [0.0], color="purple", lw=2, linestyle=":", label= L"\alpha_p = 2.4")
    ax1.plot([0.0], [0.0], color="purple", lw=2, linestyle=ls(".-"), label=L"\alpha_p = 2.2")

    # handles, labels = ax1.get_legend_handles_labels()
    # ax1.add_artist(legend(handles[1:2], labels[1:2], frameon=false, loc="lower left"))
    # ax1.add_artist(legend(handles[3:5], labels[3:5], frameon=false, loc="upper right"))

    ax1.legend(frameon=false, loc="lower left", bbox_to_anchor=(-0.9, -0.8), ncols=4)

    #legend(loc="upper right", frameon=false)

    savefig(plot_name, bbox_inches="tight")
    close(fig)
end


names = ["Coma", "Virgo", "Perseus"]
plot_name = plot_path * "Fig05.pdf"

plot_Fig05(names, plot_name)

