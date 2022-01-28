using JuMP
using CPLEX
include("parser.jl")


function plans_coupants(MyFileName::String)
  n, s, t, S, d1, d2, p, ph, d, D = read_instance(MyFileName)

  # Create the model
  m = Model(CPLEX.Optimizer)

  # probleme maitre
  ## Variables
  @variable(m, x[1:n,1:n], Bin) # certains arcs n'existent pas, ils ne sont pas utilises dans les contraintes
  @variable(m, y[1:n], Bin)
  @variable(m, z >= 0)

  ## Constraints
  @constraint(m, sum(x[s,j] for j in 1:n if d[s,j] != 0) == 1)
  @constraint(m, sum(x[i,t] for i in 1:n if d[i,t] != 0) == 1)
  @constraint(m, [v in 1:n; v != s && v != t], sum(x[i,v] for i in 1:n if d[i,v] != 0) == sum(x[v,j] for j in 1:n if d[v,j] != 0))
  @constraint(m, [v in 1:n; v != t], y[v] == sum(x[v,j] for j in 1:n if d[v,j] != 0))
  @constraint(m, [i in 1:n; i != t], y[t] == sum(x[i,t] for i in 1:n if d[i,t] != 0))
  @constraint(m, sum(p[i]*y[i] for i in 1:n) <= S)
  @constraint(m, z >= sum(d[i,j]*x[i,j] for i in 1:n, j in 1:n if d[i,j] != 0)) # objectif robuste reformule

  ## Objective
  @objective(m, Min, z)

  #resolution
  optimize!(m)

  z_star = JuMP.objective_value.(m)
  x_star = JuMP.getvalue.( m[:x] )
  y_star = JuMP.getvalue.( m[:y] )

  z1 = z_star + 1
  z2 = S + 1

  while z2 >= S + 1e-3 || z1 > z_star + 1e-3 && z1 < z_star - 1e-3

    #### sous-probleme 1 ####
    # Create the model
    m1 = Model(CPLEX.Optimizer)

    ## Variables
    @variable(m1, delta1[1:n,1:n] >= 0)

    ## Constraints
    @constraint(m1,[i in 1:n, j in 1:n], delta1[i,j] <= D[i,j])
    @constraint(m1, sum(delta1[i,j] for i in 1:n, j in 1:n) <= d1)

    ## Objective
    @objective(m1, Max, sum(d[i,j]*(1 + delta1[i,j])*x_star[i,j] for i in 1:n, j in 1:n if d[i,j] != 0))

    #resolution
    optimize!(m1)
    z1 = JuMP.objective_value.(m1)
    delta_1 = JuMP.getvalue.( m1[:delta1] )


    #### sous-probleme 2 ####
    # Create the model
    m2 = Model(CPLEX.Optimizer)

    ## Variables
    @variable(m2, delta2[1:n] >= 0)

    ## Constraints
    @constraint(m2,[v in 1:n], delta2[v] <= 2)
    @constraint(m2, sum(delta2[v] for v in 1:n) <= d2)

    ## Objective
    @objective(m2, Max, sum((p[i] + delta2[i]*ph[i])*y_star[i] for i in 1:n))

    #resolution
    optimize!(m2)
    z2 = JuMP.objective_value.(m2)
    delta_2 = JuMP.getvalue.( m2[:delta2] )

    if z2 >= S + 1e-3 || z1 > z_star + 1e-3 && z1 < z_star - 1e-3 # pour gerer la premiere iteration
      #### adding Constraints
      @constraint(m, z >= sum(x[i,j]*d[i,j]*(1+ delta_1[i,j]) for i in 1:n, j in 1:n if d[i,j] != 0))
      @constraint(m, sum(y[v]*(p[v] + ph[v]*delta_2[v]) for v in 1:n) <= S)

      optimize!(m)
      z_star = JuMP.objective_value.(m)
    end
  end
end
