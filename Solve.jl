include("Include.jl")

# load the data dictionary -
#data_dictionary = maximize_acetate_data_dictionary(0,0,0)
#data_dictionary = maximize_atp_data_dictionary(0,0,0)
data_dictionary = maximize_cellmass_data_dictionary(0,0,0)

# Turn on or turn off regulation by commenting out this line
#include("Regulation.jl")
# solve the lp problem -
(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)

list = show_flux_profile_markdown(flux_array,0.00001,data_dictionary)
writedlm("anaerobic.md",list)
println(flux_array[24])
