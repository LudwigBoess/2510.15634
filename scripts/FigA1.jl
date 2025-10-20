include(joinpath(@__DIR__, "config.jl"))

using GadgetUnits
using SpectralCRsUtility
using PyPlot, PyPlotUtility
using DelimitedFiles

"""
    Spectrum
"""

global const x_inj = 3.3
global const mp = 1.6726e-24
global const cL = 2.9979e10
global const kB = 1.38066e-16
global const GeVtoerg = 1.60217733e-3
global const slope_soft = 1.e-6



"""
    find_init_norm(pressure::T, slope::T, bound_low::T, bound_up::T) where T

Find norm of first bin for a given total pressure.
"""
function find_init_norm(energy_density::T, slope::T, bound_low::T, bound_up::T) where T

    norm_cnst = energy_density / ( 4π * cL * bound_low^4 )

    init_norm = norm_cnst * ( 4 - slope ) / ( (bound_up/bound_low)^(4 - slope) - 1 )

    if ( 4 - slope_soft ) < slope < ( 4 + slope_soft )

        slope_var = (slope - 4)/slope_soft
        init_norm2 = norm_cnst / log(bound_up/bound_low)

        if !iszero(slope_var)
            return init_norm * slope_var + init_norm2 * ( 1 - slope_var )
        else
            return init_norm2
        end
    end

    return init_norm
end


function plot_FigA1(Emin, plot_name)

    PyPlot.matplotlib[:rc]("text", usetex=true)
    GU = GadgetPhysical()
    nH = 1.0
    rho = nH / GU.rho_ncm3
    rho_cgs = rho * GU.rho_cgs

    ucr = 1.0 * GeVtoerg / GU.P_CR_cgs
    U = 1.0 / GU.T_eV
    T = 1.e8 # K

    m_p = 1.6726e-24
    m_e = 9.10953e-28
    k_B = 1.38066e-16
    c_light = 2.9979e10

    Nbins = 100

    p0 = Emin * GeVtoerg / (mp * cL^2)
    pmax = 1.e5 * GeVtoerg / (mp * cL^2)
    par = CRMomentumDistributionConfig(p0, pmax, 10)
    bounds = momentum_bin_boundaries(par)
    
    Eγ = 10.0 .^ LinRange(-4, 2, Nbins)
    jγ = Vector{Float64}(undef, Nbins)

    cmap = PyPlot.cm.plasma
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.Normalize(vmin=1, vmax=5))
    sm.set_array([])

    colors = [color=sm.to_rgba(i) for i = 1:4]
    q0 = [2.2, 2.4, 2.7, 3.0]

    x_pixels = 700
    fig = get_figure(1.0; x_pixels)
    plot_styling!(x_pixels)

    gs = plt.GridSpec(2, 1, figure=fig, height_ratios=[1, 0.3], hspace=0.1)

    subplot(get_gs(gs, 0, 0))

    ax1 = gca()
    axis_ticks_styling!(ax1)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim([1.e-4, 1.e2])
    #ax1.set_ylim([1.e-22, 3.e-15])
    ax1.set_ylim([1.e-22, 8.e-15])
    ax1.set_ylabel( L"q_\gamma(E_\gamma)" * " [" * L"\gamma" * " cm" * L"^{-3}" * " s" * L"^{-1}" * " GeV" * L"^{-1}" * "]")
    ax1.set_xticklabels([])


    subplot(get_gs(gs, 1, 0))

    ax2 = gca()
    axis_ticks_styling!(ax2)
    ax2.set_xscale("log")

    ax2.set_xlim([1.e-4, 1.e2])
    ax2.set_ylim([-0.5, 0.5])
    ax2.set_ylabel(L"\Delta q_\gamma/q_\gamma (E_\gamma)")
    ax2.set_xlabel(L"\gamma" * "-photon Energy " * L"E_\gamma" * " [GeV]")

    for i = 1:length(q0)

        jγ_minot = readdlm(joinpath(data_path, "appendix/spectrum_$(Emin)GeV_$(q0[i]).txt"))
        #jγ_minot = readdlm(joinpath(@__DIR__, "spectrum_10GeV_$(q0[i]).txt"))

        # plot analytic solution
        ax1.plot(Eγ, jγ_minot, color=colors[i], linestyle="--")

        q = q0[i] + 2
        f_0 = find_init_norm(ucr, q, p0, pmax)
        fp = [f_0 * 10.0^(j * par.bin_width * (-q)) for j = 0:par.Nbins-1] .* GU.CR_norm 

        for j = 1:length(Eγ)
            jγ[j] = gamma_source_pions(fp, q .* ones(par.Nbins), pmax, bounds, nH, Eγ[j],
                                        N_subcycle = 10)
        end

        println("error: ", (maximum(jγ) / maximum(jγ_minot)))

        ax1.plot(Eγ, jγ, color=colors[i], label=L"q = %$(q)")

        ax2.plot(Eγ, (jγ .- jγ_minot) ./ jγ_minot, color=colors[i])
    end

    ax1.plot([0.0], [0.0], label=L"\texttt{minot}", color="k", linestyle="--")
    ax1.plot([0.0], [0.0], label=L"\texttt{Crescendo}", color="k")
    ax1.legend(frameon=false)

    savefig(plot_name, bbox_inches="tight", dpi=300)
    close(fig)
end

Emin = 1
plot_name = plot_path * "FigA1a.pdf"
plot_FigA1(Emin, plot_name)

Emin = 10
plot_name = plot_path * "FigA1b.pdf"
plot_FigA1(Emin, plot_name)

