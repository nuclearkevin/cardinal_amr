#--------------------------------------------------------------------------------------------------------------------------#
# This is based on the C5G7 reactor physics benchmark problem extension case as described in:                              #
# "Benchmark on Deterministic Transport Calculations Without Spatial Homogenisation: MOX Fuel Assembly 3-D Extension Case" #
# [NEA/NSC/DOC(2005)16]                                                                                                    #
# https://www.oecd-nea.org/upload/docs/application/pdf/2019-12/nsc-doc2005-16.pdf                                          #
#                                                                                                                          #
# The original C5G7 benchmark is defined with multi-group cross sections. To account for                                   #
# continuous energy spectral effects, we chose to use the material properties provided in:                                 #
# "Proposal for a Second Stage of the Benchmark on Power Distributions Within Assemblies"                                  #
# [NEA/NSC/DOC(96)2]                                                                                                       #
# https://www.oecd-nea.org/upload/docs/application/pdf/2020-01/nsc-doc96-02-rev2.pdf                                       #
#--------------------------------------------------------------------------------------------------------------------------#

import openmc
import numpy as np
import openmc_common as geom
from openmc_materials import MATERIALS as mats

#--------------------------------------------------------------------------------------------------------------------------#
# Geometry definitions.
## All of the pincells.
PINCELLS = {}

## The pincells - fuels first
### Common primitives for defining the different fuel regions.
fuel_pin_or = openmc.ZCylinder(r = geom.r_fuel)
fuel_gap_or = openmc.ZCylinder(r = geom.r_fuel + geom.t_f_c_gap)
fuel_zr_or = openmc.ZCylinder(r = geom.r_fuel + geom.t_f_c_gap + geom.t_zr_clad)
fuel_bb = openmc.model.RectangularPrism(width = geom.pitch, height = geom.pitch)

### The entire 4.3% MOX pincell.
gap_cell_1 = openmc.Cell(name = '4.3% MOX Pin Gap 1', region = +fuel_pin_or & -fuel_gap_or)
zr_clad_cell_1 = openmc.Cell(name = '4.3% MOX Pin Zr Clad', region = +fuel_gap_or & -fuel_zr_or, fill = mats['ZR_C'])
h2o_bb_cell_1 = openmc.Cell(name = '4.3% MOX Pin Water Bounding Box', region = +fuel_zr_or & -fuel_bb, fill = mats['H2O'])
mox_4_3_fuel_cell = openmc.Cell(name = '4.3% MOX Fuel Pin', region = -fuel_pin_or, fill = mats['MOX_43'])
PINCELLS['MOX43'] = openmc.Universe(cells=[mox_4_3_fuel_cell, gap_cell_1, zr_clad_cell_1, h2o_bb_cell_1])

### The entire 7.0% MOX pincell.
gap_cell_2 = openmc.Cell(name = '7.0% MOX Pin Gap 1', region = +fuel_pin_or & -fuel_gap_or)
zr_clad_cell_2 = openmc.Cell(name = '7.0% MOX Pin Zr Clad', region = +fuel_gap_or & -fuel_zr_or, fill = mats['ZR_C'])
h2o_bb_cell_2 = openmc.Cell(name = '7.0% MOX Pin Water Bounding Box', region = +fuel_zr_or & -fuel_bb, fill = mats['H2O'])
mox_7_0_fuel_cell = openmc.Cell(name = '7.0% MOX Fuel Pin', region = -fuel_pin_or, fill = mats['MOX_70'])
PINCELLS['MOX70'] = openmc.Universe(cells=[mox_7_0_fuel_cell, gap_cell_2, zr_clad_cell_2, h2o_bb_cell_2])

### The entire 8.7% MOX pincell.
gap_cell_3 = openmc.Cell(name = '8.7% MOX Pin Gap 1', region = +fuel_pin_or & -fuel_gap_or)
zr_clad_cell_3 = openmc.Cell(name = '8.7% MOX Pin Zr Clad', region = +fuel_gap_or & -fuel_zr_or, fill = mats['ZR_C'])
h2o_bb_cell_3 = openmc.Cell(name = '8.7% MOX Pin Water Bounding Box', region = +fuel_zr_or & -fuel_bb, fill = mats['H2O'])
mox_8_7_fuel_cell = openmc.Cell(name = '8.7% MOX Fuel Pin', region = -fuel_pin_or, fill = mats['MOX_87'])
PINCELLS['MOX87'] = openmc.Universe(cells=[mox_8_7_fuel_cell, gap_cell_3, zr_clad_cell_3, h2o_bb_cell_3])

### The entire UO2 pincell.
gap_cell_4 = openmc.Cell(name = 'UO2 Pin Gap 1', region = +fuel_pin_or & -fuel_gap_or)
zr_clad_cell_4 = openmc.Cell(name = 'UO2 Pin Zr Clad', region = +fuel_gap_or & -fuel_zr_or, fill = mats['ZR_C'])
h2o_bb_cell_4 = openmc.Cell(name = 'UO2 Pin Water Bounding Box', region = +fuel_zr_or & -fuel_bb, fill = mats['H2O'])
uo2_fuel_cell = openmc.Cell(name = 'UO2 Fuel Pin', region = -fuel_pin_or, fill = mats['UO2'])
PINCELLS['UO2'] = openmc.Universe(cells=[uo2_fuel_cell, gap_cell_4, zr_clad_cell_4, h2o_bb_cell_4])

## Guide tube, control rod, and fission chamber next.
### Common primitives for defining both.
tube_fill_or = openmc.ZCylinder(r = geom.r_guide)
tube_clad_or = openmc.ZCylinder(r = geom.r_guide + geom.t_al_clad)

tube_clad_cell_1 = openmc.Cell(name = 'Guide Tube Cladding', region = +tube_fill_or & -tube_clad_or, fill = mats['AL_C'])
tube_clad_cell_2 = openmc.Cell(name = 'Control Rod Cladding', region = +tube_fill_or & -tube_clad_or, fill = mats['AL_C'])
tube_clad_cell_3 = openmc.Cell(name = 'Fission Chamber Cladding', region = +tube_fill_or & -tube_clad_or, fill = mats['AL_C'])

guide_tube_h2o_bb_cell_1 = openmc.Cell(name = 'Guide Tube Water Bounding Box', region = +tube_clad_or & -fuel_bb, fill = mats['H2O'])
guide_tube_h2o_bb_cell_2 = openmc.Cell(name = 'Control Rod Water Bounding Box', region = +tube_clad_or & -fuel_bb, fill = mats['H2O'])
guide_tube_h2o_bb_cell_3 = openmc.Cell(name = 'Fission Chamber Water Bounding Box', region = +tube_clad_or & -fuel_bb, fill = mats['H2O'])

### The guide tube.
tube_fill_cell = openmc.Cell(name = 'Guide Tube Water', region = -tube_fill_or, fill = mats['H2O'])
PINCELLS['GTB'] = openmc.Universe(cells=[tube_fill_cell, tube_clad_cell_1, guide_tube_h2o_bb_cell_1])

### The control rod.
rod_fill_cell = openmc.Cell(name = 'Control Rod Meat', region = -tube_fill_or, fill = mats['BC4'])
PINCELLS['ROD'] = openmc.Universe(cells=[rod_fill_cell, tube_clad_cell_2, guide_tube_h2o_bb_cell_2])

### The fission chamber.
fission_chamber_cell = openmc.Cell(name = 'Fission Chamber', region = -tube_fill_or, fill = mats['FISS'])
PINCELLS['FIS'] = openmc.Universe(cells=[fission_chamber_cell, tube_clad_cell_3, guide_tube_h2o_bb_cell_3])

### An empty water block.
water_cell = openmc.Cell(name = 'Moderator', region = -fuel_bb, fill = mats['H2O'])
PINCELLS['H2O'] = openmc.Universe(cells=[water_cell])
#--------------------------------------------------------------------------------------------------------------------------#PINC
