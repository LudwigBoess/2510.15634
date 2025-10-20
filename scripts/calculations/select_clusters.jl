using GadgetIO, GadgetUnits
using DelimitedFiles

"""
    filter_masses(filename)

Filter function to use with GadgetIO.jl. 
Searches a subfind file for all clusters with Mvir > 5.e14.
"""
function filter_masses(filename)

    GU = GadgetPhysical(read_header(filename))
    mvir = read_subfind(filename, "M200") .* GU.m_msun

    if length(mvir) > 0
        return findall(mvir .> 5.0e13)
    else
        return Int64[]
    end
end

sub_base = "/e/ocean3/Local/3072/nonrad_mhd_crs_new/groups_000_z=0/sub_000"
clusters_M_gt_5e13 = filter_subfind(sub_base, filter_masses)

function write_halo_ids(halo_ids, filename)
    
    data = Matrix{Int64}(undef, 2, length(halo_ids))
    for i in 1:length(halo_ids)
        data[1, i] = halo_ids[i].file
        data[2, i] = halo_ids[i].id
    end

    writedlm(filename, data')
end

filename = "/e/ocean2/users/lboess/PaperRepos/GammaWeb/data/clusters_M200_gt_5e13_new.txt"
write_halo_ids(clusters_M_gt_5e13, filename)

#println(clusters_M_gt_5e13)