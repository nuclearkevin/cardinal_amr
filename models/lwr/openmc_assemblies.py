#--------------------------------------------------------------------------------------------------------------------------#
# This is a modified version of the C5G7 reactor physics benchmark problem extension case as described in:                 #
# "Benchmark on Deterministic Transport Calculations Without Spatial Homogenisation: MOX Fuel Assembly 3-D Extension Case" #
# [NEA/NSC/DOC(2005)16]                                                                                                    #
# https://www.oecd-nea.org/upload/docs/application/pdf/2019-12/nsc-doc2005-16.pdf                                          #                                                                                                                     #
#                                                                                                                          #
# The original C5G7 benchmark is defined with multi-group cross sections. To account for                                   #
# continuous energy spectral effects, we chose to use the material properties provided in:                                 #
# "Proposal for a Second Stage of the Benchmark on Power Distributions Within Assemblies"                                  #
# [NEA/NSC/DOC(96)2]                                                                                                       #
# https://www.oecd-nea.org/upload/docs/application/pdf/2020-01/nsc-doc96-02-rev2.pdf                                       #
#--------------------------------------------------------------------------------------------------------------------------#

import openmc
import openmc_common as geom
import openmc_pincells as p

#--------------------------------------------------------------------------------------------------------------------------#
# Geometry definitions.
## The assemblies.
assembly_bb = openmc.model.RectangularPrism(width = 17.0 * geom.pitch, height = 17.0 * geom.pitch)

### UO2 fueled assembly.
uo2_assembly = openmc.RectLattice(name = 'UO2 Assembly')
uo2_assembly.pitch = (geom.pitch, geom.pitch)
uo2_assembly.lower_left = (-17.0 * geom.pitch / 2.0, -17.0 * geom.pitch / 2.0)
uo2_assembly.universes = [
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 1
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 2
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 3
  [p.uo2_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.uo2_u], # 4
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 5
  [p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u], # 6
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 7
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 8
  [p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.fis_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u], # 9
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 10
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 11
  [p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u], # 12
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 13
  [p.uo2_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.uo2_u], # 14
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.tub_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 15
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 16
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u]  # 17
]# 1        2        3        4        5        6        7        8        9        10       11       12       13       14       15       16       17
uo2_assembly_uni = openmc.Universe(cells = [openmc.Cell(name = 'UO2 Assembly Cell', region = -assembly_bb, fill = uo2_assembly)])

### UO2 fueled assembly with inserted control rods.
uo2_rodded_assembly = openmc.RectLattice(name = 'Rodded UO2 Assembly')
uo2_rodded_assembly.pitch = (geom.pitch, geom.pitch)
uo2_rodded_assembly.lower_left = (-17.0 * geom.pitch / 2.0, -17.0 * geom.pitch / 2.0)
uo2_rodded_assembly.universes = [
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 1
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 2
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 3
  [p.uo2_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.uo2_u], # 4
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 5
  [p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u], # 6
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 7
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 8
  [p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.fis_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u], # 9
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 10
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 11
  [p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u], # 12
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 13
  [p.uo2_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.uo2_u], # 14
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.rod_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 15
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u], # 16
  [p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u, p.uo2_u]  # 17
]# 1        2        3        4        5        6        7        8        9        10       11       12       13       14       15       16       17
uo2_rodded_assembly_uni = openmc.Universe(cells = [openmc.Cell(name = 'Rodded UO2 Assembly Cell', region = -assembly_bb, fill = uo2_rodded_assembly)])

### MOX fueled assembly.
mox_assembly = openmc.RectLattice(name = 'MOX Assembly')
mox_assembly.pitch = (geom.pitch, geom.pitch)
mox_assembly.lower_left = (-17.0 * geom.pitch / 2.0, -17.0 * geom.pitch / 2.0)
mox_assembly.universes = [
  [p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u], # 1
  [p.mox43_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox43_u], # 2
  [p.mox43_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.tub_u,   p.mox70_u, p.mox70_u, p.tub_u,   p.mox70_u, p.mox70_u, p.tub_u,   p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox43_u], # 3
  [p.mox43_u, p.mox70_u, p.mox70_u, p.tub_u,   p.mox70_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox70_u, p.tub_u,   p.mox70_u, p.mox70_u, p.mox43_u], # 4
  [p.mox43_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox43_u], # 5
  [p.mox43_u, p.mox70_u, p.tub_u,   p.mox87_u, p.mox87_u, p.tub_u,   p.mox87_u, p.mox87_u, p.tub_u,   p.mox87_u, p.mox87_u, p.tub_u,   p.mox87_u, p.mox87_u, p.tub_u,   p.mox70_u, p.mox43_u], # 6
  [p.mox43_u, p.mox70_u, p.mox70_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox70_u, p.mox70_u, p.mox43_u], # 7
  [p.mox43_u, p.mox70_u, p.mox70_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox70_u, p.mox70_u, p.mox43_u], # 8
  [p.mox43_u, p.mox70_u, p.tub_u,   p.mox87_u, p.mox87_u, p.tub_u,   p.mox87_u, p.mox87_u, p.fis_u,   p.mox87_u, p.mox87_u, p.tub_u,   p.mox87_u, p.mox87_u, p.tub_u,   p.mox70_u, p.mox43_u], # 9
  [p.mox43_u, p.mox70_u, p.mox70_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox70_u, p.mox70_u, p.mox43_u], # 10
  [p.mox43_u, p.mox70_u, p.mox70_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox70_u, p.mox70_u, p.mox43_u], # 11
  [p.mox43_u, p.mox70_u, p.tub_u,   p.mox87_u, p.mox87_u, p.tub_u,   p.mox87_u, p.mox87_u, p.tub_u,   p.mox87_u, p.mox87_u, p.tub_u,   p.mox87_u, p.mox87_u, p.tub_u,   p.mox70_u, p.mox43_u], # 12
  [p.mox43_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox43_u], # 13
  [p.mox43_u, p.mox70_u, p.mox70_u, p.tub_u,   p.mox70_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox87_u, p.mox70_u, p.tub_u,   p.mox70_u, p.mox70_u, p.mox43_u], # 14
  [p.mox43_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.tub_u,   p.mox70_u, p.mox70_u, p.tub_u,   p.mox70_u, p.mox70_u, p.tub_u,   p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox43_u], # 15
  [p.mox43_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox70_u, p.mox43_u], # 16
  [p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u, p.mox43_u]  # 17
]# 1          2          3          4          5          6          7          8          9          10         11         12         13         14         15         16         17
mox_assembly_uni = openmc.Universe(cells = [openmc.Cell(name = 'MOX Assembly Cell', region = -assembly_bb, fill = mox_assembly)])

### The portion of the upper reflector containing control rods.
rodded_ref_assembly = openmc.RectLattice(name = 'Rodded Reflector Assembly')
rodded_ref_assembly.pitch = (geom.pitch, geom.pitch)
rodded_ref_assembly.lower_left = (-17.0 * geom.pitch / 2.0, -17.0 * geom.pitch / 2.0)
rodded_ref_assembly.universes = [
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 1
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 2
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 3
  [p.h2o_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.h2o_u], # 4
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 5
  [p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u], # 6
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 7
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 8
  [p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u], # 9
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 10
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 11
  [p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u], # 12
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 13
  [p.h2o_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.h2o_u], # 14
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.rod_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 15
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 16
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u]  # 17
]# 1        2        3        4        5        6        7        8        9        10       11       12       13       14       15       16       17
rodded_ref_assembly_uni = openmc.Universe(cells = [openmc.Cell(name = 'Rodded Reflector Assembly Cell', region = -assembly_bb, fill = rodded_ref_assembly)])

### The portion of the upper reflector containing guide tubes.
unrodded_ref_assembly = openmc.RectLattice(name = 'Unrodded Reflector Assembly')
unrodded_ref_assembly.pitch = (geom.pitch, geom.pitch)
unrodded_ref_assembly.lower_left = (-17.0 * geom.pitch / 2.0, -17.0 * geom.pitch / 2.0)
unrodded_ref_assembly.universes = [
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 1
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 2
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 3
  [p.h2o_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.h2o_u], # 4
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 5
  [p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u], # 6
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 7
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 8
  [p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u], # 9
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 10
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 11
  [p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u], # 12
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 13
  [p.h2o_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.h2o_u], # 14
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.tub_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 15
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u], # 16
  [p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u, p.h2o_u]  # 17
]# 1        2        3        4        5        6        7        8        9        10       11       12       13       14       15       16       17
unrodded_ref_assembly_uni = openmc.Universe(cells = [openmc.Cell(name = 'Unrodded Reflector Assembly Cell', region = -assembly_bb, fill = unrodded_ref_assembly)])
#--------------------------------------------------------------------------------------------------------------------------#
