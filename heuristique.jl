using JuMP
using CPLEX
include("parser.jl")

function heuristique(MyFileName::String)
  n, s, t, S, d1, d2, p, ph, d, D = read_instance(MyFileName)

  # Create the model
  m = Model(CPLEX.Optimizer)

end
