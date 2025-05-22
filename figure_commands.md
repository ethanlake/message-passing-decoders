## commands for reproducing plots

### offline decoding

**1D**
plog vs p: 
`simulation_plotter.py -fin sim_data/1dloc_erode_L32_smallrand.jld2 sim_data/1dloc_erode_L64_smallrand.jld2 sim_data/1dloc_erode_L128_smallrand.jld2 sim_data/1dloc_erode_L256_smallrand.jld2 sim_data/1dloc_erode_L512_smallrand.jld2 -just_f`

plog vs L:
`simulation_plotter.py -fin sim_data/1dloc_erode_p0.3_fscaling.jld2 sim_data/1dloc_erode_p0.35_fscaling.jld2 sim_data/1dloc_erode_p0.4_fscaling.jld2 -just_f`

decoding time: 
`simulation_plotter.py -fin sim_data/1dloc_erode_p0.2.jld2 sim_data/1dloc_erode_p0.25.jld2 sim_data/1dloc_erode_p0.3.jld2 -just_te`

uncoordinated plog vs p: 
`simulation_plotter.py -fin sim_data/1dloc_erode_L32_asynch.jld2 sim_data/1dloc_erode_L64_asynch.jld2 sim_data/1dloc_erode_L128_asynch.jld2 sim_data/1dloc_erode_L256_asynch.jld2 sim_data/1dloc_erode_L512_asynch.jld2 -just_f`

**2D** 
plog vs p: 
`simulation_plotter.py -fin sim_data/2dloc_erode_L8.jld2 sim_data/2dloc_erode_L16.jld2 sim_data/2dloc_erode_L32.jld2 sim_data/2dloc_erode_L64.jld2 sim_data/2dloc_erode_L128.jld2  -just_f`

plog vs L: 
`simulation_plotter.py -fin sim_data/2dloc_erode_p0.045_fscaling.jld2 sim_data/2dloc_erode_p0.05_fscaling.jld2 sim_data/2dloc_erode_p0.055_fscaling.jld2 -just_f`

decoding time: 
`simulation_plotter.py -fin sim_data/2dloc_erode_p0.01.jld2 sim_data/2dloc_erode_p0.02.jld2 sim_data/2dloc_erode_p0.03.jld2 -just_te`

uncoordinated plog vs p: 
`simulation_plotter.py -fin sim_data/2dloc_erode_L8_asynch.jld2 sim_data/2dloc_erode_L16_asynch.jld2 sim_data/2dloc_erode_L32_asynch.jld2 sim_data/2dloc_erode_L64_asynch.jld2 sim_data/2dloc_erode_L128_asynch.jld2 -just_f`

uncoordinated plog vs L: 
`simulation_plotter.py -fin sim_data/2dloc_erode_p0.03_asynch_fscaling.jld2  sim_data/2dloc_erode_p0.035_asynch_fscaling.jld2 sim_data/2dloc_erode_p0.04_asynch_fscaling.jld2 -just_f`

### online decoding 

**windowed decoding**

**direct message passing**

1D tmem vs L for direct msg passing: 
`simulation_plotter.py -fin sim_data/1dloc_trel_p0.15.jld2 sim_data/1dloc_trel_p0.135.jld2 sim_data/1dloc_trel_p0.125.jld2`

1D tmem vs p for direct msg passing: 
` simulation_plotter.py -fin sim_data/1dloc_trel_L8.jld2 sim_data/1dloc_trel_L16.jld2 sim_data/1dloc_trel_L32.jld2 sim_data/1dloc_trel_L64.jld2 sim_data/1dloc_trel_L128.jld2`

**instantaneous power-law interactions**

alpha = 1/4: 
`simulation_plotter.py -fin sim_data/powerdw_trel_alpha0.25_L16_gam0.jld2 sim_data/powerdw_trel_alpha0.25_L32_gam0.jld2 sim_data/powerdw_trel_alpha0.25_L64_gam0.jld2 sim_data/powerdw_trel_alpha0.25_L128_gam0.jld2`

alpha = 3/4:
`simulation_plotter.py -fin sim_data/powerdw_trel_alpha0.75_L8_gam0.jld2 sim_data/powerdw_trel_alpha0.75_L16_gam0.jld2 sim_data/powerdw_trel_alpha0.75_L32_gam0.jld2 sim_data/powerdw_trel_alpha0.75_L64_gam0.jld2 sim_data/powerdw_trel_alpha0.75_L128_gam0.jld2`

alpha = 2: 
`simulation_plotter.py -fin sim_data/powerdw_trel_alpha2.0_L8_gam0.jld2 sim_data/powerdw_trel_alpha2.0_L16_gam0.jld2 sim_data/powerdw_trel_alpha2.0_L32_gam0.jld2 sim_data/powerdw_trel_alpha2.0_L64_gam0.jld2 sim_data/powerdw_trel_alpha2.0_L128_gam0.jld2`
