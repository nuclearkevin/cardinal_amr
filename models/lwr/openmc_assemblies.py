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
from openmc_pincells import PINCELLS as pins

#--------------------------------------------------------------------------------------------------------------------------#
# Geometry definitions.
## The assemblies.
pins_per_axis = 17.0
ASSEMBLIES = {}
assembly_bb = openmc.model.RectangularPrism(width = pins_per_axis * geom.pitch, height = pins_per_axis * geom.pitch)

### Different assembly maps.
uo2_assembly_map = [
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 1
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 2
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 3
  ['UO2', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'UO2'], # 4
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 5
  ['UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2'], # 6
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 7
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 8
  ['UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'FIS', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2'], # 9
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 10
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 11
  ['UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2'], # 12
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 13
  ['UO2', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'UO2'], # 14
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'GTB', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 15
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 16
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2']  # 17
]#  1      2      3      4      5      6      7      8      9      10     11     12     13     14     15     16     17

rodded_uo2_assembly_map = [
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 1
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 2
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 3
  ['UO2', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'UO2'], # 4
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 5
  ['UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2'], # 6
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 7
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 8
  ['UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'FIS', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2'], # 9
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 10
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 11
  ['UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2'], # 12
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 13
  ['UO2', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'UO2'], # 14
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'ROD', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 15
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2'], # 16
  ['UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2', 'UO2']  # 17
]#  1      2      3      4      5      6      7      8      9      10     11     12     13     14     15     16     17

mox_assembly_map = [
  ['MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43'], # 1
  ['MOX43', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX43'], # 2
  ['MOX43', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'GTB',   'MOX70', 'MOX70', 'GTB',   'MOX70', 'MOX70', 'GTB',   'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX43'], # 3
  ['MOX43', 'MOX70', 'MOX70', 'GTB',   'MOX70', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX70', 'GTB',   'MOX70', 'MOX70', 'MOX43'], # 4
  ['MOX43', 'MOX70', 'MOX70', 'MOX70', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX70', 'MOX70', 'MOX70', 'MOX43'], # 5
  ['MOX43', 'MOX70', 'GTB',   'MOX87', 'MOX87', 'GTB',   'MOX87', 'MOX87', 'GTB',   'MOX87', 'MOX87', 'GTB',   'MOX87', 'MOX87', 'GTB',   'MOX70', 'MOX43'], # 6
  ['MOX43', 'MOX70', 'MOX70', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX70', 'MOX70', 'MOX43'], # 7
  ['MOX43', 'MOX70', 'MOX70', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX70', 'MOX70', 'MOX43'], # 8
  ['MOX43', 'MOX70', 'GTB',   'MOX87', 'MOX87', 'GTB',   'MOX87', 'MOX87', 'FIS',   'MOX87', 'MOX87', 'GTB',   'MOX87', 'MOX87', 'GTB',   'MOX70', 'MOX43'], # 9
  ['MOX43', 'MOX70', 'MOX70', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX70', 'MOX70', 'MOX43'], # 10
  ['MOX43', 'MOX70', 'MOX70', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX70', 'MOX70', 'MOX43'], # 11
  ['MOX43', 'MOX70', 'GTB',   'MOX87', 'MOX87', 'GTB',   'MOX87', 'MOX87', 'GTB',   'MOX87', 'MOX87', 'GTB',   'MOX87', 'MOX87', 'GTB',   'MOX70', 'MOX43'], # 12
  ['MOX43', 'MOX70', 'MOX70', 'MOX70', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX70', 'MOX70', 'MOX70', 'MOX43'], # 13
  ['MOX43', 'MOX70', 'MOX70', 'GTB',   'MOX70', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX87', 'MOX70', 'GTB',   'MOX70', 'MOX70', 'MOX43'], # 14
  ['MOX43', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'GTB',   'MOX70', 'MOX70', 'GTB',   'MOX70', 'MOX70', 'GTB',   'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX43'], # 15
  ['MOX43', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX70', 'MOX43'], # 16
  ['MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43', 'MOX43']  # 17
]#  1        2        3        4        5        6        7        8        9        10       11       12       13       14       15       16       17

rodded_refl_assembly_map = [
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 1
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 2
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 3
  ['H2O', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'H2O'], # 4
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 5
  ['H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O'], # 6
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 7
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 8
  ['H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O'], # 9
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 10
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 11
  ['H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O'], # 12
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 13
  ['H2O', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'H2O'], # 14
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'ROD', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 15
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 16
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O']  # 17
]#  1      2      3      4      5      6      7      8      9      10     11     12     13     14     15     16     17

refl_assembly_map = [
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 1
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 2
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 3
  ['H2O', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'H2O'], # 4
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 5
  ['H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O'], # 6
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 7
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 8
  ['H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O'], # 9
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 10
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 11
  ['H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O'], # 12
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 13
  ['H2O', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'H2O'], # 14
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'GTB', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 15
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O'], # 16
  ['H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O']  # 17
]#  1      2      3      4      5      6      7      8      9      10     11     12     13     14     15     16     17

### UO2 fueled assembly.
uo2_assembly = openmc.RectLattice(name = 'UO2 Assembly')
uo2_assembly.pitch = (geom.pitch, geom.pitch)
uo2_assembly.lower_left = (-pins_per_axis * geom.pitch / 2.0, -pins_per_axis * geom.pitch / 2.0)
uo2_assembly.universes = [ [ pins[pin_type] for pin_type in row] for row in uo2_assembly_map ]
ASSEMBLIES['UO2'] = openmc.Universe(cells = [openmc.Cell(name = 'UO2 Assembly Cell', region = -assembly_bb, fill = uo2_assembly)])

### UO2 fueled assembly with inserted control rods.
uo2_rodded_assembly = openmc.RectLattice(name = 'Rodded UO2 Assembly')
uo2_rodded_assembly.pitch = (geom.pitch, geom.pitch)
uo2_rodded_assembly.lower_left = (-pins_per_axis * geom.pitch / 2.0, -pins_per_axis * geom.pitch / 2.0)
uo2_rodded_assembly.universes = [ [ pins[pin_type] for pin_type in row] for row in rodded_uo2_assembly_map ]
ASSEMBLIES['UO2_ROD'] = openmc.Universe(cells = [openmc.Cell(name = 'Rodded UO2 Assembly Cell', region = -assembly_bb, fill = uo2_rodded_assembly)])

### MOX fueled assembly.
mox_assembly = openmc.RectLattice(name = 'MOX Assembly')
mox_assembly.pitch = (geom.pitch, geom.pitch)
mox_assembly.lower_left = (-pins_per_axis * geom.pitch / 2.0, -pins_per_axis * geom.pitch / 2.0)
mox_assembly.universes = [ [ pins[pin_type] for pin_type in row] for row in mox_assembly_map ]
ASSEMBLIES['MOX'] = openmc.Universe(cells = [openmc.Cell(name = 'MOX Assembly Cell', region = -assembly_bb, fill = mox_assembly)])

### The portion of the upper reflector containing control rods.
rodded_ref_assembly = openmc.RectLattice(name = 'Rodded Reflector Assembly')
rodded_ref_assembly.pitch = (geom.pitch, geom.pitch)
rodded_ref_assembly.lower_left = (-pins_per_axis * geom.pitch / 2.0, -pins_per_axis * geom.pitch / 2.0)
rodded_ref_assembly.universes = [ [ pins[pin_type] for pin_type in row] for row in rodded_refl_assembly_map ]
ASSEMBLIES['REF_ROD'] = openmc.Universe(cells = [openmc.Cell(name = 'Rodded Reflector Assembly Cell', region = -assembly_bb, fill = rodded_ref_assembly)])

### The portion of the upper reflector containing guide tubes.
unrodded_ref_assembly = openmc.RectLattice(name = 'Unrodded Reflector Assembly')
unrodded_ref_assembly.pitch = (geom.pitch, geom.pitch)
unrodded_ref_assembly.lower_left = (-pins_per_axis * geom.pitch / 2.0, -pins_per_axis * geom.pitch / 2.0)
unrodded_ref_assembly.universes = [ [ pins[pin_type] for pin_type in row] for row in refl_assembly_map ]
ASSEMBLIES['REF'] = openmc.Universe(cells = [openmc.Cell(name = 'Unrodded Reflector Assembly Cell', region = -assembly_bb, fill = unrodded_ref_assembly)])
#--------------------------------------------------------------------------------------------------------------------------#
