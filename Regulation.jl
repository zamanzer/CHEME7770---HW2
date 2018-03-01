

flux = data_dictionary["default_flux_bounds_array"]

species = data_dictionary["species_bounds_array"]

(size_boolean,) = size(data_dictionary["default_flux_bounds_array"])
control = trues(size_boolean)

regulators = trues(21)
index = 1
while index < 11
########### Regulators #######################
# ArcA - 1
if species[89,1] == 0
regulators[1] = 1
else
regulators[1] = 0
end

# DcuR - 2
if regulators[3] != 0
regulators[2] = 1
else
regulators[2] = 0
end

# DcuS - 3
if species[92,1] != 0 || species[80,1] != 0 || species[87,1] != 0
regulators[3] = 1
else
regulators[3] = 0
end

# FadR - 4
if species[81,1] != 0 || species[73,1] == 0
regulators[4] = 1
else
regulators[4] = 0
end

# Fis - 5
#always on

# FnR - 6
if species[89,1] == 0
regulators[6] = 1
else
regulators[6] = 0
end

# FruR - 7
if regulators[20] == 0
regulators[7] = 1
else
regulators[7] = 0
end

# GlcC - 8
if species[73,1] != 0
regulators[8] = 1
else
regulators[8] = 0
end

# GlnG - 9
if species[88,1] == 0
regulators[9] = 1
else
regulators[9] = 0
end

# IcIR - 10
if regulators[4] != 0
regulators[10] = 1
else
regulators[10] = 0
end

# Mlc - 11
if species[81,1] == 0
regulators[11] = 1
else
regulators[11] = 0
end

# Nac - 12
if regulators[19] != 0
regulators[12] = 1
else
regulators[12] = 0
end

# PdhR - 13
if regulators[21] == 0
regulators[13] = 1
else
regulators[13] = 0
end

# PhoB - 14
if regulators[15] != 0
regulators[14] = 1
else
regulators[14] = 0
end

# PhoR - 15
if species[90,1] == 0
regulators[15] = 1
else
regulators[15] = 0
end

# CRPnoGLC - 16
if species[81,1] == 0
regulators[16] = 1
else
regulators[16] = 0
end

# CRPnoGLM - 17
temp1 = species[81,1] + species[87,1] + species[86,1]
if  temp1 == 0
regulators[17] = 1
else
regulators[17] = 0
end

# NRI_hi - 18
if regulators[19] != 0
regulators[18] = 1
else
regulators[18] = 0
end

# NRI_low - 19
if regulators[9] != 0
regulators[19] = 1
else
regulators[19] = 0
end

# surplusFDP - 20
#temp2 = flux[139,2] + flux[140,2] + flux[134,2] + flux[135,2] + flux[108,2] + flux[109,2] #139,134,108 have reverse
#if (flux[64,2] == 0 && temp2 == 0) || species[79,1] != 0
#regulators[20] = 1
#else
#regulators[20] = 0
#end

# surplusPYR - 21
temp3 = flux[98,2] + flux[97,2]
temp4 = flux[76,2] + flux[122,2] + flux[106,2] + flux[91,2] + flux[92,2] + flux[129,2] #91 has reverse
if temp3 == 0 && temp4 == 0
regulators[21] = 1
else
regulators[21] = 0
end

##################################################

#################FLUX CONTROLLER##################

# ICL - 90
temp5 = min(1-regulators[1],regulators[7])
if regulators[10] == 0 && temp5 == 0
control[90] = 1
else
control[90] = 0
end

# MALS - 93
if control[90] != 0
control[93] = 1
else
control[93] = 0
end

# PDH - 105
if regulators[13] == 0 || regulators[5] != 0
control[105] = 1
else
control[105] = 0
end

# ACALD - 1(has reverse) and ALCD2x - 18 (has reverse)
temp6 = min(species[89,1],regulators[7])
if max(regulators[5],species[89,1],1-temp6) != 0
control[1] = 1
control[2] = 1
control[18] = 1
control[19] = 1
else
control[1] = 0
control[2] = 0
control[18] = 0
control[19] = 0
end

# CYTBD - 28
if regulators[6] == 0  || regulators[1] != 0
control[28] = 1
else
control[28] = 0
end

# FUMt2_2 - 71 and MALt2_2 - 94 and SUCCt2_2 - 129
if regulators[17] != 0 && regulators[1] == 0 && regulators[2] != 0
control[71] = 1
control[94] = 1
control[129] = 1
else
control[71] = 0
control[94] = 0
control[129] = 0
end

# FORt2 - 65 and FORti - 66
temp7 = regulators[6] + regulators[16]
if regulators[1] != 0 || temp7 == 2
control[65] = 1
control[66] = 1
else
control[65] = 0
control[66] = 0
end

# FRD7 - 67
if regulators[6] != 0 || regulators[2] != 0
control[67] = 1
else
control[67] = 0

end

# FUM - 69 has reverse
temp8 = max(regulators[1],regulators[6])
if 1-temp8 == 0
control[69] = 1
control[70] = 1
else
control[69] = 0
control[70] = 0
end

if regulators[6] != 0 || regulators[16] != 0 || regulators[2] != 0
control[69] = 1
control[70] = 1
else
control[69] = 0
control[70] = 0
end

if regulators[1] == 0
control[69] = 1
control[70] = 1
else
control[69] = 0
control[70] = 0
end

# GLUDy - 79 (has reverse)
if regulators[12] == 0 && species[83,1] == 0
control[79] = 1
control[80] = 1
else
control[79] = 0
control[80] = 0
end

# D_LACt2 - 29 (has reverse)
if regulators[1] == 0 && regulators[8] != 0
control[29] = 1
control[30] = 1
else
control[29] = 0
control[30] = 0
end

# MALS - 93
if control[29] != 0
control[93] = 1
else
control[93] = 0
end

# GLNS - 77 (Changed per Prof. Varner recommendation)
if regulators[16] == 0
control[77] = 1
else
control[77] = 0
end

# GLUSy - 82
if regulators[18] == 0 || species[83,1] == 0
control[82] = 1
else
control[82] = 0
end

# D_LACt2 - 29 (has reverse)
if regulators[1] == 0
control[29] = 1
control[30] = 1
else
control[29] = 0
control[30] = 0
end

# FRUpts2 - 68 and GLCpts - 76
if regulators[17] != 0 || regulators[11] == 0
control[68] = 1
control[76] = 1
else
control[68] = 0
control[76] = 0
end

# MDH - 95 (has reverse)
if regulators[1] == 0
control[95] = 1
control[96] = 1
else
control[95] = 0
control[96] = 0
end

# NADH16 - 99
if regulators[1] == 0 && regulators[6] == 0
control[99] = 1
else
control[99] = 0
end

# PFL - 107
temp9 = min(regulators[6],regulators[16])
if regulators[1] != 0 || temp9 != 0
control[107] = 1
else
control[107] = 0
end

if regulators[1] != 0 || regulators[6] != 0
control[107] = 1
else
control[107] = 0
end

# PIt2r - 115 (has reverse)
if regulators[14] == 0
control[115] = 1
control[116] = 1
else
control[115] = 0
control[116] = 0
end

# PPS - 119
if regulators[7] != 0
control[119] = 1
else
control[119] = 0
end

# GLCpts - 76
if regulators[11] == 0 || regulators[7] == 0
control[76] = 1
else
control[76] = 0
end

# PYK - 122
if regulators[7] == 0
control[122] = 1
else
control[122] = 0
end

# SUCDi - 131
temp10 = max(regulators[1],regulators[6])
if regulators[5] != 0 || regulators[16] != 0 || 1-temp10 != 0
control[131] = 1
else
control[131] = 0
end

# ACKr - 5 (has reverse)
if regulators[16] != 0 || regulators[6] != 0
control[5] = 1
control[6] = 1
else
control[5] = 0
control[6] = 0
end

# PFL - 107
if control[5] != 0
control[107] = 1
else
control[107] = 0
end

# GLUN - 81
temp10 = min(species[55,1],1-regulators[16])
if species[35,1] == 0 || temp10 != 0
control[81] = 1
else
control[81] = 0
end
index = index + 1

end

flux[:,2] = flux[:,2].*control
data_dictionary["default_flux_bounds_array"] = flux
###########################################################
