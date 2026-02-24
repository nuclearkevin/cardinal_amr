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

import sys
sys.path.append("../")

import numpy as np
from argparse import ArgumentParser

import openmc
import openmc_common as geom
from openmc_materials import MATERIALS as mats
import openmc_pincells as pins
from openmc_settings import COMMON_SETTINGS as settings

ap = ArgumentParser()
ap.add_argument('-n', dest='n_axial', type=int, default=50,
                help='Number of axial cell divisions')
args = ap.parse_args()

#--------------------------------------------------------------------------------------------------------------------------#
# Geometry definitions.
pins.fuel_bb.boundary_type = 'reflective'
fuel_bot = openmc.ZPlane(z0 = -geom.core_height)
fuel_top = openmc.ZPlane(z0 = geom.core_height)

core_assembly = openmc.RectLattice(name = 'UO2 Pincell Lattice')
core_assembly.pitch = (geom.pitch, geom.pitch, 2.0 * geom.core_height / args.n_axial)
core_assembly.lower_left = (-geom.pitch / 2.0, -geom.pitch / 2.0, -geom.core_height)
core_assembly.universes = [ [ [pins.PINCELLS['UO2'] ] ] for i in range(args.n_axial) ]
core_assembly_cell = openmc.Cell(name = 'UO2 Pincell Lattice Cell',
                                 region = -pins.fuel_bb & +fuel_bot & -fuel_top,
                                 fill = core_assembly)

## Add the top and bottom axial water reflector.
## Set the boundary condition on the topmost plane to vacuum.
refl_top = openmc.ZPlane(z0 = geom.core_height + geom.reflector_t, boundary_type = 'vacuum')
refl_top_cell = openmc.Cell(name='Axial Reflector Cell (Top)', fill = mats['H2O'], region=-pins.fuel_bb & -refl_top & +fuel_top)
refl_bot = openmc.ZPlane(z0 = -geom.core_height - geom.reflector_t, boundary_type = 'vacuum')
refl_bot_cell = openmc.Cell(name='Axial Reflector Cell (Bottom)', fill = mats['H2O'], region=-pins.fuel_bb & +refl_bot & -fuel_bot)
#--------------------------------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------------------------------#
# Setup the model.
pincell_model = openmc.Model(geometry = openmc.Geometry(openmc.Universe(cells = [core_assembly_cell, refl_top_cell, refl_bot_cell])))

## The simulation settings.
pincell_model.settings = settings
pincell_model.settings.source = [openmc.IndependentSource(space = openmc.stats.Box(lower_left = (-geom.pitch / 2.0, -geom.pitch / 2.0, -geom.core_height),
                                                                                   upper_right = (geom.pitch / 2.0, geom.pitch / 2.0, geom.core_height)))]

pincell_model.export_to_model_xml()
#--------------------------------------------------------------------------------------------------------------------------#
