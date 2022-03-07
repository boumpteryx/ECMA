include("parser.jl")
include("dualisation.jl")
include("branch_and_cut.jl")
include("heuristique.jl")
include("plans_coupants.jl")

start = time()

### resolution par dualisation
dualisation("20_USA-road-d.BAY.gr")

### resolution par plans coupants
# plans_coupants("40_USA-road-d.COL.gr")

### resolution par LazyCallback (branch-and-cut)
# branch_and_cut("20_USA-road-d.COL.gr")

### heuristique
# heuristique("20_USA-road-d.BAY.gr")

files_and_dirs = readdir("./temporaire")
p = Vector{Char}(undef,0)
secs = Vector{Float64}(undef,0)
for file in files_and_dirs
    actual_time = time()
    println("filename = ", file)
    result = plans_coupants(file)
    println("filename = ", file)
    append!(p, file)
    local computation_time = time() - actual_time
    append!(secs, computation_time)
    println(computation_time, " seconds")
end
println(p)
println(secs)
