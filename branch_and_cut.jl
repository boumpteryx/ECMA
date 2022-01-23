using JuMP
using CPLEX
include("parser.jl")



function branch_and_cut(MyFileName::String)
  n, s, t, S, d1, d2, p, ph, d, D = read_instance(MyFileName)

  # Create the model
  m = Model(CPLEX.Optimizer)

  # Il est imposé d'utiliser 1 seul thread en Julia avec CPLEX pour utiliser les callbacks
  MOI.set(m, MOI.NumberOfThreads(), 1)

end
