include("parser.jl")
include("dualisation.jl")
include("branch_and_cut.jl")
include("heuristique.jl")
include("plans_coupants.jl")

### resolution par dualisation
dualisation("testyMcTestFace")

### resolution par plans coupants
# plans_coupants("testyMcTestFace")

### resolution par LazyCallback (branch-and-cut)
# branch_and_cut("testyMcTestFace")

### heuristique
# heuristique("testyMcTestFace")
