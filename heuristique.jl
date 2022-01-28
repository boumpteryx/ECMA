using JuMP
using CPLEX
include("parser.jl")

function dual_relax(MyFileName::String)
  n, s, t, S, d1, d2, p, ph, d, D = read_instance(MyFileName)

  # Create the model
  m = Model(CPLEX.Optimizer)

  ## Variables
  # x relaxe, y relaxe
  @variable(m, x[1:n,1:n] >= 0) # certains arcs n'existent pas, ils ne sont pas utilises dans les contraintes
  @variable(m, y[1:n] >= 0)
  @variable(m, t1 >= 0)
  @variable(m, t2 >= 0)
  @variable(m, z[1:n,1:n] >= 0, Int) # meme remarque
  @variable(m, zprim[1:n] >= 0, Int)

  ## Constraints
  @constraint(m,[i in 1:n], y[i] <= 1) # relaxation
  @constraint(m,[i in 1:n, j in 1:n], x[i,j] <= 1) # relaxation
  @constraint(m,sum(x[s,j] for j in 1:n if d[s,j] != 0) == 1)
  @constraint(m,sum(x[i,t] for i in 1:n if d[i,t] != 0) == 1)
  @constraint(m,[v in 1:n; v != s && v != t], sum(x[i,v] for i in 1:n if d[i,v] != 0) == sum(x[v,j] for j in 1:n if d[v,j] != 0))
  @constraint(m,[v in 1:n; v != t], y[v] == sum(x[v,j] for j in 1:n if d[v,j] != 0))
  @constraint(m,[i in 1:n; i != t], y[t] == sum(x[i,t] for i in 1:n if d[i,t] != 0))
  @constraint(m, t2*d2 + sum(p[i]*y[i] + 2*zprim[i] for i in 1:n) <= S)
  @constraint(m, [i in 1:n, j in 1:n ; d[i,j] != 0], t1 + z[i,j] >= d[i,j]*x[i,j])
  @constraint(m, [i in 1:n], t2 + zprim[i] >= ph[i]*y[i])

  ## Objective
  @objective(m, Min, d1*t1 + sum(d[i,j]*x[i,j] + D[i,j]*z[i,j] for i in 1:n, j in 1:n if d[i,j] != 0))

  #resolution
  optimize!(m)

  return JuMP.objective_value.(m) # enlever . et JuMP
end

function heuristic(MyFileName::String)
  n, s, t, S, d1, d2, p, ph, d, D = read_instance(MyFileName)

  relax_dual_value = dual_relax(MyFileName)

  # do heuristic


end
