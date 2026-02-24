!include ../common.i

T_fluid = ${fparse 280.0 + 273.15}

[Mesh]
  [file]
    type = FileMeshGenerator
    file = mesh_in.e
  []
[]

[Variables]
  [temp]
    initial_condition = ${T_fluid}
  []
[]

[AuxVariables]
  [heat_source]
    family = MONOMIAL
    order = CONSTANT
    block = '0 1'
  []
  [heat_source_rel_error]
    family = MONOMIAL
    order = CONSTANT
    block = '0 1'
  []
[]

[Adaptivity]
  marker = error_frac
  steps = 1

  [Indicators/error]
    type = GradientJumpIndicator
    variable = temp
    block = '0 1'
  []
  [Markers/error_frac]
    type = ErrorFractionLookAheadMarker
    # Statistical error
    rel_error_refine = 0.01
    stat_error_indicator = 'heat_source_rel_error'
    # Spatial error
    indicator = error
    refine = 0.3
    coarsen = 0.0
  []
[]

[Kernels]
  [hc]
    type = HeatConduction
    variable = temp
    block = '0 1 2 3'
  []
  [heat]
    type = CoupledForce
    variable = temp
    v = heat_source
    block = '0 1'
  []
[]

[Functions]
  [T_fluid]
    type = ParsedFunction
    expression = '573.0 + 50.0 * ((z + ${core_height}) / ${fparse 2.0 * core_height})'
  []
[]

[BCs]
  [surface]
    type = ConvectiveFluxFunction
    T_infinity = T_fluid

    # convert from W/m2/K to W/cm2/K
    coefficient = ${fparse 1000.0/100.0/100.0}
    variable = temp
    boundary = 'clad_or'
  []
[]

[Materials]
  [k_clad]
    type = GenericConstantMaterial
    prop_values = '0.5'
    prop_names = 'thermal_conductivity'
    block = '3'
  []
  [k_fuel]
    type = GenericConstantMaterial
    prop_values = '0.05'
    prop_names = 'thermal_conductivity'
    block = '0 1'
  []
  [k_gap]
    type = GenericConstantMaterial
    prop_values = '0.1'
    prop_names = 'thermal_conductivity'
    block = '2'
  []
[]

[Executioner]
  type = Transient
  nl_abs_tol = 1e-8
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
  csv = true
[]

[Postprocessors]
  [source_integral]
    type = ElementIntegralVariablePostprocessor
    variable = heat_source
    execute_on = transfer
    block = '0 1'
  []
  [max_T]
    type = NodalExtremeValue
    variable = temp
  []
[]
