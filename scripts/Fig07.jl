"""
    Explanatory figure for emissivity scaling with slope and p_inj
"""

using SpectralCRsUtility
using GadgetUnits
using Base.Threads
using PyPlot, PyPlotUtility
using ProgressMeter
using DelimitedFiles

global const slope_soft = 1.e-4
global const x_inj = 3.3
global const mp = 1.6726e-24
global const cL = 2.9979e10
global const kB = 1.38066e-16
global const GeVtoerg = 1.60217733e-3


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


function get_phase_map(p_range, q_range, Nbins)

    GU = GadgetPhysical()

    B = 1.0e-6
    pmax = 1.e6
    #ϵ_cr = 1.0
    #ϵ_cr_GU = ϵ_cr * GU.x_cgs^3 / GU.E_cgs
    ϵ_cr_GU = 1.0 * GeVtoerg / GU.P_CR_cgs


    j_γ = Matrix{Float64}(undef, Nbins, Nbins)

    @showprogress for i = 1:Nbins
        for j = 1:Nbins

            p_min = p_range[i]
            q = q_range[j]

            par = CRMomentumDistributionConfig(p_min, pmax, 24)
            bounds = momentum_bin_boundaries(par)

            f_0 = find_init_norm(ϵ_cr_GU, q, p_min, pmax)
            fp = [f_0 * 10.0^(j * par.bin_width * (-q)) for j = 0:par.Nbins-1] .* GU.CR_norm

            fq = q .* ones(par.Nbins)
            cut = pmax

            j_γ[i, j] = gamma_luminosity_pions(fp, fq, cut, bounds, 1.0, 1.0)
        end
    end

    return j_γ
end

function plot_Fig07(q_range, p_range, j_nu, plot_name)

    c_lim = [1.e-23, 2.e-19]
    levels = 10.0 .^ LinRange(log10(c_lim[1]), log10(c_lim[2]), 
        floor(Int64,abs((log10(c_lim[2]) - log10(c_lim[1])) * 2)))

    fig = get_figure(1.0)
    plot_styling!()
    ax = gca()
    ax.set_yscale("log")
    xlabel("Spectral Slope  " * L"q")
    ylabel("Injection Momentum  " * L"\hat{p}_\mathrm{inj}" * "  [" * L"m_p c" * "]")
    im = pcolormesh(q_range, p_range, j_nu,
        cmap="plasma", norm=matplotlib.colors.LogNorm(vmin=c_lim[1], vmax=c_lim[2]))

    CS = contour(q_range, p_range, j_nu, levels=levels, colors="white", 
        linewidths=1.0, alpha=0.5)
    CS = contour(q_range, p_range, j_nu, levels=levels, colors="k", 
        linewidths=0.5, linestyles="dashed")
    #ax.clabel(CS, CS.levels, inline=true)

    get_colorbar_right(ax, im, "γ-ray Emissivity  " * L"j_\gamma" * " [erg s" * L"^{-1}" * "cm" * L"^{-3}" * "]")

    savefig(plot_name, bbox_inches="tight", dpi=500)
    close(fig)
end

Nbins = 100
p_range = 10.0 .^ LinRange(-2, 0, Nbins)
q_range = LinRange(4, 6, Nbins)

filename = data_path * "gamma_data/emissivity.txt"

#j_γ = get_phase_map(p_range, q_range, Nbins)
#writedlm(filename, j_γ)

j_γ = readdlm(filename)

println(maximum(j_γ))


plot_name = plot_path * "Fig07.pdf"

plot_Fig07(q_range, p_range, j_γ, plot_name)



