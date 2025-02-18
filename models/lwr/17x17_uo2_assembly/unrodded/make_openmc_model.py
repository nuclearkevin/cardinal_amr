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
sys.path.append("../../")

import numpy as np
from argparse import ArgumentParser

import openmc
import openmc_common as geom
import openmc_materials as mats
import openmc_assemblies as asmb

ap = ArgumentParser()
ap.add_argument('-n', dest='n_axial', type=int, default=1,
                help='Number of axial core divisions')
args = ap.parse_args()

#--------------------------------------------------------------------------------------------------------------------------#
# Geometry definitions.
core_z_planes = [ openmc.ZPlane(z0=z) for z in np.linspace(0.0, geom.core_height, args.n_axial + 1) ]
core_z_planes[0].boundary_type = 'reflective'

all_cells = []
for layer_idx, planes in enumerate(zip(core_z_planes[:-1], core_z_planes[1:])):
  all_cells.append(openmc.Cell(name = f'UO2 Assembly Cell {layer_idx}', region = -asmb.assembly_bb & +planes[0] & -planes[1], fill = asmb.uo2_assembly))

refl_top = openmc.ZPlane(z0 = geom.core_height + geom.reflector_t, boundary_type = 'vacuum')
all_cells.append(openmc.Cell(name='Axial Reflector Cell', fill = asmb.unrodded_ref_assembly_uni, region=-asmb.assembly_bb & -refl_top & +core_z_planes[-1]))

#--------------------------------------------------------------------------------------------------------------------------#
# Setup the model.
c5g7_model = openmc.Model(geometry = openmc.Geometry(openmc.Universe(cells = all_cells)), materials = openmc.Materials([mats.uo2, mats.h2o, mats.fiss, mats.zr, mats.al]))

## The simulation settings.
c5g7_model.settings.source = [openmc.IndependentSource(space = openmc.stats.Box(lower_left = (-17.0 * geom.pitch / 2.0, -17.0 * geom.pitch / 2.0, 0.0),
                                                                                upper_right = (17.0 * geom.pitch / 2.0, 17.0 * geom.pitch / 2.0, geom.core_height)))]
c5g7_model.settings.batches = 100
c5g7_model.settings.generations_per_batch = 10
c5g7_model.settings.inactive = 10
c5g7_model.settings.particles = 1000

c5g7_model.export_to_model_xml()
#--------------------------------------------------------------------------------------------------------------------------#
