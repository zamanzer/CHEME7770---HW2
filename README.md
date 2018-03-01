### Generated files

This code is run by using the command include("Solve.jl") in the julia terminal.

The upper bound on the exhance species were all set to 100 in order to recreate the growth rates given in the Palsson paper. While these species don't seem logical, it agrees with the Palsson paper.
To turn on anaerobic growth calculations, either comment or uncomment the line in data_dictionary_cell_growth that turns the lower oxygen bound to 0.
To turn on or off regulations either comment or uncomment the command include("Regulation.jl") in the Solve.jl script
The values used for the actual calculations for the metabolic constraints were obtained from the following:
http://kirschner.med.harvard.edu/files/bionumbers/ComparisonoffermentationproductsofE.coliandK.oxytoca.pdf

Using these parameters for the upper bound in the anaerobic case and then setting them to 0 for the aerobic will likely give more realistic results.

As can be seen from the flux distributions given in the files "aerobic", "aerobic_regulation", etc. the growth rate does not change with the presence of gene regulation on the metabolic pathway. This signifies that the cell is largely good at maximizing cell growth on its own.

The regulation file calculates the boolean rules given in the Palsson paper and then applies them to the flux of a specific reaction.

For the Aerobic case the following fluxes/reactions are inactivated by the control, when using the appropriate bounding species: FORti, FORt2, ICL, FUMt2_2, MALS, SUCCT2_2, PFL, and MALt2_2.

For the Anaerobic case the same reactions are shut off as well as MDH, NADH16, and D_LaACt2.

