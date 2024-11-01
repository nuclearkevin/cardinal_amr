[Mesh]
  [file]
    type = FileMeshGenerator
    file = mesh_in.e
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  particles = 20000
  inactive_batches = 500
  batches = 10000

  verbose = true
  power = ${fparse 3000e6 / 273 / (17 * 17) * 9}
  cell_level = 1
  normalize_by_global_tally = false

  [Tallies]
    [heat_source]
      type = MeshTally
      score = 'kappa_fission flux'
      name = heat_source
      output = 'unrelaxed_tally_rel_error'
    []
  []
[]

[Executioner]
  type = Steady
[]

[Postprocessors]
  [num_active]
    type = NumElems
    elem_filter = active
  []
  [num_total]
    type = NumElems
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
[]