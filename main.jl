include("parser.jl")
include("dualisation.jl")


### resolution par dualisation
dualisation("testyMcTestFace")

### resolution par plans coupants

### resolution par LazyCallback
# Il est imposé d'utiliser 1 seul thread en Julia avec CPLEX pour utiliser les callbacks
MOI.set(m, MOI.NumberOfThreads(), 1)
