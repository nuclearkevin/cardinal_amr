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
from argparse import ArgumentParser

ap = ArgumentParser()
ap.add_argument('-n', dest='n_axial', type=int, default=1,
                help='Number of axial cell divisions')
args = ap.parse_args()

#--------------------------------------------------------------------------------------------------------------------------#
# Some geometric properties that can be modified to change the model.
## The radius of a fuel pin (same for all pin types).
r_fuel = 0.4095

## The thickness of the fuel-clad gap.
t_f_c_gap = 0.0085

## The thickness of the Zr fuel pin cladding.
t_zr_clad = 0.057

## The pitch of a single lattice element.
pitch = 1.26

## The height of the fuel assemblies from the axial midplane.
core_height = 192.78

## The thickness of the axial reflector above the lattice.
reflector_t = 21.42

# Some discretization parameters.
pin_axial_slices = args.n_axial
#--------------------------------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------------------------------#
# Material definitions.
## Fuel region: UO2 at ~1% enriched.
uo2_comp = 1.0e24 * np.array([8.65e-4, 2.225e-2, 4.622e-2])
uo2_frac = uo2_comp / np.sum(uo2_comp)
uo2 = openmc.Material(name = 'UO2 Fuel', temperature = 293.15)
uo2.add_nuclide('U235', uo2_frac[0], percent_type = 'ao')
uo2.add_nuclide('U238', uo2_frac[1], percent_type = 'ao')
uo2.add_element('O', uo2_frac[2], percent_type = 'ao')
uo2.set_density('atom/cm3', np.sum(uo2_comp))

## Moderator and coolant, boronated water.
h2o_comp = 1.0e24 * np.array([3.35e-2, 2.78e-5])
h2o_frac = h2o_comp / np.sum(h2o_comp)
h2o = openmc.Material(name = 'H2O Moderator', temperature = 293.15)
h2o.add_element('H', 2.0 * h2o_frac[0], percent_type = 'ao')
h2o.add_element('O', h2o_frac[0], percent_type = 'ao')
h2o.add_element('B', h2o_frac[1], percent_type = 'ao')
h2o.set_density('atom/cm3', np.sum(h2o_comp))
h2o.add_s_alpha_beta('c_H_in_H2O')

## Zr clad.
zr = openmc.Material(name = 'Zr Cladding', temperature = 293.15)
zr.add_element('Zr', 1.0, percent_type = 'ao')
zr.set_density('atom/cm3', 1.0e24 * 4.30e-2)
#--------------------------------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------------------------------#
# Geometry definitions.
## Create cylindrical surfaces
fuel_or = openmc.ZCylinder(r = r_fuel, name = 'Fuel Outer Radius')
clad_ir = openmc.ZCylinder(r = r_fuel + t_f_c_gap, name = 'Clad Inner Radius')
clad_or = openmc.ZCylinder(r = r_fuel + t_f_c_gap + t_zr_clad, name = 'Clad Outer Radius')

## Create a region represented as the inside of a rectangular prism.
box = openmc.model.RectangularPrism(pitch, pitch, boundary_type='reflective')

## Create cells, mapping materials to regions - split up the axial height.
plane_surfaces = [ openmc.ZPlane(z0=z) for z in np.linspace(0.0, core_height, pin_axial_slices + 1) ]

## set the boundary condition on the bottommost plane to reflective.
plane_surfaces[0].boundary_type = 'reflective'

## Create the actual pincells.
all_cells = []
for layer_idx, planes in enumerate(zip(plane_surfaces[:-1], plane_surfaces[1:])):
  layer = +planes[0] & -planes[1]
  all_cells.append(openmc.Cell(fill=uo2, region=-fuel_or & layer, name=f'Fuel {layer_idx}'))
  all_cells.append(openmc.Cell(fill=None, region=+fuel_or & -clad_ir & layer, name=f'Gap {layer_idx}'))
  all_cells.append(openmc.Cell(fill=zr, region=+clad_ir & -clad_or & layer, name=f'Clad {layer_idx}'))
  all_cells.append(openmc.Cell(fill=h2o, region=+clad_or & layer & -box, name=f'Water {layer_idx}'))

## Add the top axial water reflector.
## Set the boundary condition on the topmost plane to vacuum.
refl_top = openmc.ZPlane(z0 = core_height + reflector_t, boundary_type = 'vacuum')
all_cells.append(openmc.Cell(name='Axial Reflector Cell', fill = h2o, region=-box & -refl_top & +plane_surfaces[-1]))
#--------------------------------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------------------------------#
# Setup the model.
pincell_model = openmc.Model(geometry = openmc.Geometry(openmc.Universe(cells = all_cells)), materials = openmc.Materials([uo2, h2o, zr]))

## The simulation settings.
pincell_model.settings.source = [openmc.IndependentSource(space = openmc.stats.Box(lower_left = (-pitch / 2.0, -pitch / 2.0, 0.0),
                                                                                   upper_right = (pitch / 2.0, pitch / 2.0, 192.78)))]

pincell_model.settings.batches = 100
pincell_model.settings.generations_per_batch = 10
pincell_model.settings.inactive = 10
pincell_model.settings.particles = 1000

pincell_model.export_to_model_xml()
#--------------------------------------------------------------------------------------------------------------------------#
