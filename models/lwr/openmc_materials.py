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

#--------------------------------------------------------------------------------------------------------------------------#
# Material definitions.
## First fuel region: 4.3% MOX.
mox_4_3_comp = 1.0e24 * np.array([5.00e-5, 2.21e-2, 1.50e-5, 5.80e-4, 2.40e-4, 9.80e-5, 5.40e-5, 1.30e-5, 4.63e-2])
mox_4_3_frac = mox_4_3_comp / np.sum(mox_4_3_comp)
mox_4_3 = openmc.Material(name = '4.3% MOX Fuel', temperature = 293.15)
mox_4_3.add_nuclide('U235', mox_4_3_frac[0], percent_type = 'ao')
mox_4_3.add_nuclide('U238', mox_4_3_frac[1], percent_type = 'ao')
mox_4_3.add_nuclide('Pu238', mox_4_3_frac[2], percent_type = 'ao')
mox_4_3.add_nuclide('Pu239', mox_4_3_frac[3], percent_type = 'ao')
mox_4_3.add_nuclide('Pu240', mox_4_3_frac[4], percent_type = 'ao')
mox_4_3.add_nuclide('Pu241', mox_4_3_frac[5], percent_type = 'ao')
mox_4_3.add_nuclide('Pu242', mox_4_3_frac[6], percent_type = 'ao')
mox_4_3.add_nuclide('Am241', mox_4_3_frac[7], percent_type = 'ao')
mox_4_3.add_element('O', mox_4_3_frac[8], percent_type = 'ao')
mox_4_3.set_density('atom/cm3', np.sum(mox_4_3_comp))

## Second fuel region: 7.0% MOX.
mox_7_0_comp = 1.0e24 * np.array([5.00e-5, 2.21e-2, 2.40e-5, 9.30e-4, 3.90e-4, 1.52e-4, 8.40e-5, 2.00e-5, 4.63e-2])
mox_7_0_frac = mox_7_0_comp / np.sum(mox_7_0_comp)
mox_7_0 = openmc.Material(name = '7.0% MOX Fuel', temperature = 293.15)
mox_7_0.add_nuclide('U235', mox_7_0_frac[0], percent_type = 'ao')
mox_7_0.add_nuclide('U238', mox_7_0_frac[1], percent_type = 'ao')
mox_7_0.add_nuclide('Pu238', mox_7_0_frac[2], percent_type = 'ao')
mox_7_0.add_nuclide('Pu239', mox_7_0_frac[3], percent_type = 'ao')
mox_7_0.add_nuclide('Pu240', mox_7_0_frac[4], percent_type = 'ao')
mox_7_0.add_nuclide('Pu241', mox_7_0_frac[5], percent_type = 'ao')
mox_7_0.add_nuclide('Pu242', mox_7_0_frac[6], percent_type = 'ao')
mox_7_0.add_nuclide('Am241', mox_7_0_frac[7], percent_type = 'ao')
mox_7_0.add_element('O', mox_7_0_frac[8], percent_type = 'ao')
mox_7_0.set_density('atom/cm3', np.sum(mox_7_0_comp))

## Third fuel region: 8.7% MOX.
mox_8_7_comp = 1.0e24 * np.array([5.00e-5, 2.21e-2, 3.00e-5, 1.16e-3, 4.90e-4, 1.90e-4, 1.05e-4, 2.50e-5, 4.63e-2])
mox_8_7_frac = mox_8_7_comp / np.sum(mox_8_7_comp)
mox_8_7 = openmc.Material(name = '8.7% MOX Fuel', temperature = 293.15)
mox_8_7.add_nuclide('U235', mox_8_7_frac[0], percent_type = 'ao')
mox_8_7.add_nuclide('U238', mox_8_7_frac[1], percent_type = 'ao')
mox_8_7.add_nuclide('Pu238', mox_8_7_frac[2], percent_type = 'ao')
mox_8_7.add_nuclide('Pu239', mox_8_7_frac[3], percent_type = 'ao')
mox_8_7.add_nuclide('Pu240', mox_8_7_frac[4], percent_type = 'ao')
mox_8_7.add_nuclide('Pu241', mox_8_7_frac[5], percent_type = 'ao')
mox_8_7.add_nuclide('Pu242', mox_8_7_frac[6], percent_type = 'ao')
mox_8_7.add_nuclide('Am241', mox_8_7_frac[7], percent_type = 'ao')
mox_8_7.add_element('O', mox_8_7_frac[8], percent_type = 'ao')
mox_8_7.set_density('atom/cm3', np.sum(mox_8_7_comp))

## Fourth fuel region: UO2 at ~1% enriched.
uo2_comp = 1.0e24 * np.array([8.65e-4, 2.225e-2, 4.622e-2])
uo2_frac = uo2_comp / np.sum(uo2_comp)
uo2 = openmc.Material(name = 'UO2 Fuel', temperature = 293.15)
uo2.add_nuclide('U235', uo2_frac[0], percent_type = 'ao')
uo2.add_nuclide('U238', uo2_frac[1], percent_type = 'ao')
uo2.add_element('O', uo2_frac[2], percent_type = 'ao')
uo2.set_density('atom/cm3', np.sum(uo2_comp))

## Control rod meat: assumed to be B-10 carbide (B4C).
bc4 = openmc.Material(name = 'Control Rod Meat', temperature = 293.15)
bc4.add_element('B', 4.0, percent_type = 'ao')
bc4.add_element('C', 1.0, percent_type = 'ao')
bc4.set_density('atom/cm3', 2.75e23)

## Moderator and coolant, boronated water.
h2o_comp = 1.0e24 * np.array([3.35e-2, 2.78e-5])
h2o_frac = h2o_comp / np.sum(h2o_comp)
h2o = openmc.Material(name = 'H2O Moderator', temperature = 293.15)
h2o.add_element('H', 2.0 * h2o_frac[0], percent_type = 'ao')
h2o.add_element('O', h2o_frac[0], percent_type = 'ao')
h2o.add_element('B', h2o_frac[1], percent_type = 'ao')
h2o.set_density('atom/cm3', np.sum(h2o_comp))
h2o.add_s_alpha_beta('c_H_in_H2O')

## Fission chamber.
fiss_comp = 1.0e24 * np.array([3.35e-2, 2.78e-5, 1.0e-8])
fiss_frac = fiss_comp / np.sum(fiss_comp)
fiss = openmc.Material(name = 'Fission Chamber', temperature = 293.15)
fiss.add_element('H', 2.0 * fiss_frac[0], percent_type = 'ao')
fiss.add_element('O', fiss_frac[0], percent_type = 'ao')
fiss.add_element('B', fiss_frac[1], percent_type = 'ao')
fiss.add_nuclide('U235', fiss_frac[2], percent_type = 'ao')
fiss.set_density('atom/cm3', np.sum(fiss_comp))
fiss.add_s_alpha_beta('c_H_in_H2O')

## Zr clad.
zr = openmc.Material(name = 'Zr Cladding', temperature = 293.15)
zr.add_element('Zr', 1.0, percent_type = 'ao')
zr.set_density('atom/cm3', 1.0e24 * 4.30e-2)

## Al clad.
al = openmc.Material(name = 'Al Cladding', temperature = 293.15)
al.add_element('Al', 1.0, percent_type = 'ao')
al.set_density('atom/cm3', 1.0e24 * 6.0e-2)
#--------------------------------------------------------------------------------------------------------------------------#
