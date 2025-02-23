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
from openmc_pincells import PINCELLS as p

#--------------------------------------------------------------------------------------------------------------------------#
# Geometry definitions.
## The assemblies.
pins_per_axis = 17.0
ASSEMBLIES = {}
assembly_bb = openmc.model.RectangularPrism(width = pins_per_axis * geom.pitch, height = pins_per_axis * geom.pitch)

### UO2 fueled assembly.
uo2_assembly = openmc.RectLattice(name = 'UO2 Assembly')
uo2_assembly.pitch = (geom.pitch, geom.pitch)
uo2_assembly.lower_left = (-pins_per_axis * geom.pitch / 2.0, -pins_per_axis * geom.pitch / 2.0)
uo2_assembly.universes = [
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 1
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 2
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 3
  [p['UO2'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['UO2']], # 4
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 5
  [p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2']], # 6
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 7
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 8
  [p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['FIS'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2']], # 9
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 10
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 11
  [p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2']], # 12
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 13
  [p['UO2'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['UO2']], # 14
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['GTB'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 15
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 16
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']]  # 17
]# 1         2         3         4         5         6         7         8         9         10        11        12        13        14        15        16        17
ASSEMBLIES['UO2'] = openmc.Universe(cells = [openmc.Cell(name = 'UO2 Assembly Cell', region = -assembly_bb, fill = uo2_assembly)])

### UO2 fueled assembly with inserted control rods.
uo2_rodded_assembly = openmc.RectLattice(name = 'Rodded UO2 Assembly')
uo2_rodded_assembly.pitch = (geom.pitch, geom.pitch)
uo2_rodded_assembly.lower_left = (-pins_per_axis * geom.pitch / 2.0, -pins_per_axis * geom.pitch / 2.0)
uo2_rodded_assembly.universes = [
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 1
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 2
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 3
  [p['UO2'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['UO2']], # 4
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 5
  [p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2']], # 6
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 7
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 8
  [p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['FIS'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2']], # 9
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 10
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 11
  [p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2']], # 12
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 13
  [p['UO2'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['UO2']], # 14
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['ROD'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 15
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']], # 16
  [p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2'], p['UO2']]  # 17
]# 1         2         3         4         5         6         7         8         9         10        11        12        13        14        15        16        17
ASSEMBLIES['UO2_ROD'] = openmc.Universe(cells = [openmc.Cell(name = 'Rodded UO2 Assembly Cell', region = -assembly_bb, fill = uo2_rodded_assembly)])

### MOX fueled assembly.
mox_assembly = openmc.RectLattice(name = 'MOX Assembly')
mox_assembly.pitch = (geom.pitch, geom.pitch)
mox_assembly.lower_left = (-pins_per_axis * geom.pitch / 2.0, -pins_per_axis * geom.pitch / 2.0)
mox_assembly.universes = [
  [p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43']], # 1
  [p['MOX43'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX43']], # 2
  [p['MOX43'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['GTB'],   p['MOX70'], p['MOX70'], p['GTB'],   p['MOX70'], p['MOX70'], p['GTB'],   p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX43']], # 3
  [p['MOX43'], p['MOX70'], p['MOX70'], p['GTB'],   p['MOX70'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX70'], p['GTB'],   p['MOX70'], p['MOX70'], p['MOX43']], # 4
  [p['MOX43'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX43']], # 5
  [p['MOX43'], p['MOX70'], p['GTB'],   p['MOX87'], p['MOX87'], p['GTB'],   p['MOX87'], p['MOX87'], p['GTB'],   p['MOX87'], p['MOX87'], p['GTB'],   p['MOX87'], p['MOX87'], p['GTB'],   p['MOX70'], p['MOX43']], # 6
  [p['MOX43'], p['MOX70'], p['MOX70'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX70'], p['MOX70'], p['MOX43']], # 7
  [p['MOX43'], p['MOX70'], p['MOX70'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX70'], p['MOX70'], p['MOX43']], # 8
  [p['MOX43'], p['MOX70'], p['GTB'],   p['MOX87'], p['MOX87'], p['GTB'],   p['MOX87'], p['MOX87'], p['FIS'],   p['MOX87'], p['MOX87'], p['GTB'],   p['MOX87'], p['MOX87'], p['GTB'],   p['MOX70'], p['MOX43']], # 9
  [p['MOX43'], p['MOX70'], p['MOX70'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX70'], p['MOX70'], p['MOX43']], # 10
  [p['MOX43'], p['MOX70'], p['MOX70'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX70'], p['MOX70'], p['MOX43']], # 11
  [p['MOX43'], p['MOX70'], p['GTB'],   p['MOX87'], p['MOX87'], p['GTB'],   p['MOX87'], p['MOX87'], p['GTB'],   p['MOX87'], p['MOX87'], p['GTB'],   p['MOX87'], p['MOX87'], p['GTB'],   p['MOX70'], p['MOX43']], # 12
  [p['MOX43'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX43']], # 13
  [p['MOX43'], p['MOX70'], p['MOX70'], p['GTB'],   p['MOX70'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX87'], p['MOX70'], p['GTB'],   p['MOX70'], p['MOX70'], p['MOX43']], # 14
  [p['MOX43'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['GTB'],   p['MOX70'], p['MOX70'], p['GTB'],   p['MOX70'], p['MOX70'], p['GTB'],   p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX43']], # 15
  [p['MOX43'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX70'], p['MOX43']], # 16
  [p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43'], p['MOX43']]  # 17
]# 1           2           3           4           5           6           7           8           9           10          11          12          13          14          15          16          17
ASSEMBLIES['MOX'] = openmc.Universe(cells = [openmc.Cell(name = 'MOX Assembly Cell', region = -assembly_bb, fill = mox_assembly)])

### The portion of the upper reflector containing control rods.
rodded_ref_assembly = openmc.RectLattice(name = 'Rodded Reflector Assembly')
rodded_ref_assembly.pitch = (geom.pitch, geom.pitch)
rodded_ref_assembly.lower_left = (-pins_per_axis * geom.pitch / 2.0, -pins_per_axis * geom.pitch / 2.0)
rodded_ref_assembly.universes = [
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 1
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 2
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 3
  [p['H2O'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['H2O']], # 4
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 5
  [p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O']], # 6
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 7
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 8
  [p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O']], # 9
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 10
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 11
  [p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O']], # 12
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 13
  [p['H2O'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['H2O']], # 14
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['ROD'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 15
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 16
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']]  # 17
]# 1         2         3         4         5         6         7         8         9         10        11        12        13        14        15        16        17
ASSEMBLIES['REF_ROD'] = openmc.Universe(cells = [openmc.Cell(name = 'Rodded Reflector Assembly Cell', region = -assembly_bb, fill = rodded_ref_assembly)])

### The portion of the upper reflector containing guide tubes.
unrodded_ref_assembly = openmc.RectLattice(name = 'Unrodded Reflector Assembly')
unrodded_ref_assembly.pitch = (geom.pitch, geom.pitch)
unrodded_ref_assembly.lower_left = (-pins_per_axis * geom.pitch / 2.0, -pins_per_axis * geom.pitch / 2.0)
unrodded_ref_assembly.universes = [
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 1
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 2
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 3
  [p['H2O'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['H2O']], # 4
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 5
  [p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O']], # 6
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 7
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 8
  [p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O']], # 9
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 10
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 11
  [p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O']], # 12
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 13
  [p['H2O'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['H2O']], # 14
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['GTB'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 15
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']], # 16
  [p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O'], p['H2O']]  # 17
]# 1         2         3         4         5         6         7         8         9         10        11        12        13        14        15        16        17
ASSEMBLIES['REF'] = openmc.Universe(cells = [openmc.Cell(name = 'Unrodded Reflector Assembly Cell', region = -assembly_bb, fill = unrodded_ref_assembly)])
#--------------------------------------------------------------------------------------------------------------------------#
