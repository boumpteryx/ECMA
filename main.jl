include("parser.jl")
include("dualisation.jl")
include("branch_and_cut.jl")
include("heuristique.jl")
include("plans_coupants.jl")

start = time()

### resolution par dualisation
dualisation("20_USA-road-d.BAY.gr")

### resolution par plans coupants
# plans_coupants("20_USA-road-d.BAY.gr")

### resolution par LazyCallback (branch-and-cut)
# branch_and_cut("20_USA-road-d.BAY.gr")

### heuristique
# heuristique("20_USA-road-d.BAY.gr")



computation_time = time() - start
println(computation_time)


# for file in readdir("./Instances_ECMA")
  # println(file)
# end
