[Mesh]
  [file]
    type = FileMeshGenerator
    file = mesh_in.e
  []
[]

[AuxVariables]
  [cell_temperature]
    family = MONOMIAL
    order = CONSTANT
  []
  [marker_var]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [cell_temperature]
    type = CellTemperatureAux
    variable = cell_temperature
  []
[]

[Adaptivity]
  marker = marker_var
  steps = 1
[]

[Problem]
  type = OpenMCCellAverageProblem
  particles = 5000
  inactive_batches = 100
  batches = 1000

  verbose = false
  power = ${fparse 3000e6 / 273 / (17 * 17)}
  source_rate_normalization = 'kappa_fission'

  cell_level = 1
  temperature_blocks = '0 1 3'
  initial_properties = 'xml'
[]

# Include common tallies.
!include ../../mesh_tallies.i
[Problem/Tallies/heat_source]
  block = '0 1'
[]

[Executioner]
  type = Transient
  num_steps = 10
[]

[MultiApps]
  [solid]
    type = TransientMultiApp
    input_files = 'solid.i'
    execute_on = timestep_end
  []
[]

[Transfers]
  [heat_source_to_solid]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = solid
    variable = heat_source
    source_variable = heat_source
    from_postprocessors_to_be_preserved = heat_source
    to_postprocessors_to_be_preserved = source_integral
  []
  [heat_rel_to_solid]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = solid
    variable = heat_source_rel_error
    source_variable = heat_source_rel_error
  []
  [temp_from_solid]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = solid
    variable = temp
    source_variable = temp
  []
  [marker_from_solid]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = solid
    variable = marker_var
    source_variable = error_frac
  []
[]

[Postprocessors]
  [heat_source]
    type = ElementIntegralVariablePostprocessor
    variable = heat_source
    execute_on = 'TIMESTEP_END TRANSFER'
    block = '0 1'
  []
  [num_active]
    type = NumElements
    elem_filter = active
  []
  [num_total]
    type = NumElements
    elem_filter = total
  []
  [max_rel_err]
    type = TallyRelativeError
    value_type = max
    tally_score = kappa_fission
  []
[]

[Outputs]
  exodus = true
  csv = true
  execute_on = 'TIMESTEP_END'
[]
