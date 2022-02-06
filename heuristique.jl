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

function Trouve_min(P,Q,d)
  min1 = -1
  min2 = -1
  dist = -1
  for p in P
    for q in Q
      if d[p,q] > 0 && (dist == -1 || d[p,q] < dist)
        min1 = p
        min2 = q
        dist = d[p,q]
      end
    end
  end
  return min1,min2,dist
end

function Dijkstra(n,s,t,d) 
  M = Array{Int64,2}(zeros(n,3))
  # 1:Distance à s
  # 2:Prédecesseur
  # 3:Nombre de prédecesseurs
  P = [s]
  Q = [i for i in 1:n if i!=s]
  for i in 1:n
    if i != s
      M[i,1] = -1
      M[i,2] = -1
    end
  end
  for q in Q
    if d[s,q] > 0
      M[q,1] = d[q,s]
      M[q,2] = s
      M[q,3] = 1
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
      if (d[q,r] > 0) && (M[q,1] + d[q,r] < M[r,1] || M[r,1] == -1)
        M[r,1] = M[q,1] + d[q,r]
        M[r,2] = q
        M[r,3] = M[q,3]+1
      end
    end
  end
  return M
end


function sortdistance(n,d)
  Sorted_d = Array{Int64}(zeros(0))
  Sorted_i = Array{Int64}(zeros(0))
  Sorted_j = Array{Int64}(zeros(0))
  for i in 1:n
    for j in 1:n
      k = 1
      dist = d[i,j]
      done = false
      while k <= size(Sorted_d)[1] && done == false
        if dist < Sorted_d[k]
          k += 1
        else
          done = true
        end
      end
      insert!(Sorted_d,k,dist)
      insert!(Sorted_i,k,i)
      insert!(Sorted_j,k,j)
    end
  end
  return Sorted_i,Sorted_j
end

function evaluatedist(L,Sorted_i,Sorted_j,d,D,d1)
  distance = 0
  for l in 1:(size(L)[1]-1)
    distance += d[L[l],L[l+1]]
  end
  remaining = d1
  k = size(L)[1] - 1
  acc = 1
  while k > 0 && remaining > 0
    i,j = Sorted_i[acc],Sorted_j[acc]
    for l in 1:(size(L)[1]-1)
      if L[l] == i && L[l+1] == j
        k -= 1
        if D[i,j] < remaining
          distance += d[i,j]*D[i,j]
          remaining -= D[i,j]
        else
          distance += d[i,j] * remaining
          remaining = 0
        end
      end
    end
    acc = acc + 1
  end
  return distance
end

function sortweight(n,ph)
  Sorted = []
  for i in 1:n
    k = 1
    done = false
    while k <= size(Sorted)[1] && done == false
      if ph[i] < ph[Sorted[k]]
        k += 1
      else
        done =true
      end
    end
    insert!(Sorted,k,i)
  end
  return Sorted
end

function evaluateweight(L,Sorted,ph,d2)
  k = size(L)[1]
  acc = 1
  remaining = d2
  weight = 0
  while k > 0 && remaining > 0
    i = Sorted[acc]
    if i in L
      k -= 1
      if 2 < remaining
        weight += 2*ph[i]
        remaining -= 2
      else
        weight += remaining * ph[i]
        remaining = 0
      end
    end
    acc += 1
  end
  return weight
end




function TotalValue(L,Value,Weight,Relativeweight,Order,Limit)
  total = 0
  weight = Array{Int64}(zeros(0))
  relativeweight = Array{Int64}(zeros(0))
  for l in L
    total += value[l]
    push!(weight,Weight[l])
    push!(relativeweight,Relativeweight[l])
  end
  total += evaluate(L,weight,relativeweight,Order,Limit)
  return total
end
  
function coeff(i,j)
  return 1
end

function heuristic(MyFileName::String)
  n, s, t, S, d1, d2, p, ph, d, D = read_instance(MyFileName)
  global_tol = 5
  instant_tol = 2
  timelimit = 60
  relax_dual_value = dual_relax(MyFileName)
#function heuristic(instant_tol,global_tol,n, s, t, S, d1, d2, p, ph, d, D,timelimit)
  Ones = Array{Int64,2}(zeros(n,n))
  for i in 1:n
    for j in 1:n
    Ones[i,j] = 1
    end
  end
  H1 = Dijkstra(n,s,t,Ones)
  H2 = Dijkstra(n,t,s,Ones)
  Length1 = Array{Int64}(zeros(n))
  for i in [1:n]
    Length1[i] = H1[i,3] + H2[i,3]
  end

  New_d = Array{Int64,2}(zeros(n,n))
  for j in 1:n
    for i in 1:n
      if d[i,j] > 0
        length = (Length1[i]+Length1[j])/2
        extra = min(floor(d1/length),coeff(i,j)*D[i,j]*d[i,j])
        New_d[i,j] = d[i,j] + extra
      end
    end
  end
  H3 = Dijkstra(n,s,t,New_d)
  H4 = Dijkstra(n,t,s,New_d)
  Length2 = Array{Int64}(zeros(n))
  for i in [1:n]
    Length2[i] = H3[i,3] + H4[i,3]
  end

  Sortedweight = sortweight(n,p)
  println(Sortedweight)
  I,J = sortdistance(n,d)
  pmin = p[Sortedweight[n]]
  pmax = p[Sortedweight[1]]
  P_lim = [S - H2[i,3]*pmin for i in 1:n]
  Distance = n * d[I[1],J[1]]
  Path = [s]
  Weight = 0

  DecisionTree = [[s,p[s],global_tol,false]]
  DecisionMatrix = [[] for i in 1:n]
  Current = s
  t0 = time()
  Started = false
  iter = 0
  while (time() - t0 < timelimit) && (Started == false || size(DecisionTree)[1] > 1)
    iter+=1
    k = size(DecisionTree)[1]
    u,tol = DecisionTree[k][1],DecisionTree[k][3]
    for v in 1:n
      if New_d[u,v] > 0 && ((H4[v,3] - H4[u,3]) <= (-1 + min(tol,instant_tol)))
        i = 1
        done = false
        while i < size(DecisionMatrix[u])[1] && done == false
          if (H4[v,3],New_d[u,v]) > (H4[DecisionMatrix[u][i],3],New_d[u,DecisionMatrix[u][i]])
            i += 1
          else 
            done = true
          end
        end
        insert!(DecisionMatrix[u],i,v)
      end
    end
    println(DecisionMatrix[s])
    forward = false
    while forward == false && (k > 1 || Started == false)
      Started = true
      while size(DecisionMatrix[u])[1] == 0 && k > 1
        pop!(DecisionTree)
        k -= 1
        if k > 0
          u,tol = DecisionTree[k][1],DecisionTree[k][3]
          if size(DecisionMatrix[u])[1] > 0
            popfirst!(DecisionMatrix[u])
          end
        else 
          Current = s
        end
      end
      if size(DecisionMatrix[u])[1] > 0
        v = DecisionMatrix[u][1]
        w = DecisionTree[k][2] + p[v]
        checkedweight = DecisionTree[k][4]
        if checkedweight == false && S - w < 2*d2*pmax
          checkedweight = true
          L = [DecisionTree[i][1] for i in 1:k]
          push!(L,v)
          w += evaluateweight(L,Sortedweight,ph,d2)
        end
        path = [DecisionTree[i][1] for i in 1:size(DecisionTree)[1]]
        if w > P_lim[v] || v in path
         popfirst!(DecisionMatrix[u])
        elseif v == t
          push!(path,v)
          weight = evaluateweight(path,Sortedweight,ph,d2)
          for i in path
            weight += p[i]
          end
          distance = evaluatedist(path,I,J,d,D,d1)
          if weight <= S && distance < Distance
            Path = path
            Weight = weight
            Distance = distance
          end
          popfirst!(DecisionMatrix[u])
        else
          push!(DecisionTree,[v,w,tol-(H4[v,3]-H4[u,3]+1),checkedweight])
          forward = true
          Current = v
        end  
      end
    end
  end
  println(Path," ",Weight," ",Distance," ",(Distance/relax_dual_value)-1," ",iter)
  return Path,Weight,Distance,(Distance/relax_dual_value)-1
end


  
  
  
  
    


