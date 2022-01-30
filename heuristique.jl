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

function Trouve_min(P,Q,d):
  min1 = -1
  min2 = -1
  dist = -1
  for p in P
    for q in Q
      if d[p][q] > 0 && (dist == -1 || d[p][q] < dist)
        min1 = q
        min2 = p
        dist = d[p][q]
      end
    end
  end
  return min1,min2,dist
end

function Dijkstra(n,s,d) 
  M = Array{Int64,3}(zeros(n,n,n))
  P = Array{Int64}(zeros(1))
  Q = Array{Int64}(zeros(n-1))
  shift = false
  for i in [1:n]
    if shift == false
      if i != s
        Q[i] = i
      else 
        P[1] = s
        shift = true
      end
    else
      Q[i-1] = i
    end
  end
  for i in [1:n]
    if i != s
      M[i,1] = -1
      M[i,2] = -1
    end
  end
  for q in Q
    if d[s][q] > 0
      M[q][1] = d[q][r]
      M[q][2] = q
      M[q][3] = 1
    end
  end
  while size(Q)[1] > 0
    p,q,dist = Trouve_min(P,Q,d)
    push!(P,q)
    j = 1
    while Q[j] != q
      j = j+1
    end
    deleteat!(Q,j)
    for r in Q
      if (d[q][r] > 0) && (M[q][1] + d[q][r] < M[r][1] || M[r][1] == -1)
        M[r][1] = M[q][1] + d[q][r]
        M[r][2] = q
        M[r][3] = M[q][3]+1
      end
    end
  end
  return M
end


function heuristic(MyFileName::String)
  n, s, t, S, d1, d2, p, ph, d, D = read_instance(MyFileName)

  relax_dual_value = dual_relax(MyFileName)

  M = Dijkstra(n,s,d)


end
