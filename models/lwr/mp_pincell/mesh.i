!include ../pincells.i

[Mesh]
  [Pin_3D]
    type = AdvancedExtruderGenerator
    input = 'UO2_Pin'
    heights = '${fparse 2.0 * core_height}'
    num_layers = '${AXIAL_DIVISIONS}'
    direction = '0.0 0.0 1.0'
  []
  [To_Origin]
    type = TransformGenerator
    input = 'Pin_3D'
    transform = TRANSLATE_CENTER_ORIGIN
  []
  [Label_Fuel_Outer]
    type = SideSetsBetweenSubdomainsGenerator
    input = 'To_Origin'
    primary_block = '1'
    paired_block = '2'
    new_boundary = 'fuel_or'
  []
  [Label_Clad_Inner]
    type = SideSetsBetweenSubdomainsGenerator
    input = 'Label_Fuel_Outer'
    primary_block = '3'
    paired_block = '2'
    new_boundary = 'clad_ir'
  []
  [Label_Clad_Outer]
    type = SideSetsBetweenSubdomainsGenerator
    input = 'Label_Clad_Inner'
    primary_block = '3'
    paired_block = '17'
    new_boundary = 'clad_or'
  []
  [Delete_Gap]
    type = BlockDeletionGenerator
    input = Label_Clad_Outer
    block = '17'
    delete_exteriors = false
  []

  final_generator = Delete_Gap
[]
