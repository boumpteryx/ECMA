using JuMP
using CPLEX
include("parser.jl")

"""
Fonction permettant de tester si un callback est appelé car une solution entière a été obtenue
"""
function isIntegerPoint(cb_data::CPLEX.CallbackContext, context_id::Clong)

    # context_id == CPX_CALLBACKCONTEXT_CANDIDATE si le callback est appelé dans un des deux cas suivants :
    # cas 1 - une solution entière a été obtenue; ou
    # cas 2 - une relaxation non bornée a été obtenue
    if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
        return false
    end

    # Pour déterminer si on est dans le cas 1 ou 2, on essaie de récupérer la
    # solution entière courante
    ispoint_p = Ref{Cint}()
    ret = CPXcallbackcandidateispoint(cb_data, ispoint_p)

    # S'il n'y a pas de solution entière
    if ret != 0 || ispoint_p[] == 0
        return false
    else
        return true
    end
end


function branch_and_cut(MyFileName::String)
  n, s, t, S, d1, d2, p, ph, d, D = read_instance(MyFileName)

  # Create the model
  m = Model(CPLEX.Optimizer)

  # Il est imposé d'utiliser 1 seul thread en Julia avec CPLEX pour utiliser les callbacks
  MOI.set(m, MOI.NumberOfThreads(), 1)

  # probleme maitre
  ## Variables
  @variable(m, x[1:n,1:n], Bin) # certains arcs n'existent pas, ils ne sont pas utilises dans les contraintes
  @variable(m, y[1:n], Bin)
  @variable(m, z >= 0)

  ## Constraints
  @constraint(m,sum(x[s,j] for j in 1:n if d[s,j] != 0) == 1)
  @constraint(m,sum(x[i,t] for i in 1:n if d[i,t] != 0) == 1)
  @constraint(m,[v in 1:n; v != s && v != t], sum(x[i,v] for i in 1:n if d[i,v] != 0) == sum(x[v,j] for j in 1:n if d[v,j] != 0))
  @constraint(m,[v in 1:n; v != t], y[v] == sum(x[v,j] for j in 1:n if d[v,j] != 0))
  @constraint(m,[i in 1:n; i != t], y[t] == sum(x[i,t] for j in 1:n if d[j,t] != 0))
  @constraint(m, sum(p[i]*y[i] for i in 1:n) <= S)
  @constraint(m, z >= sum(d[i][j]*x[i][j] for i in 1:n, j in 1:n if d[i][j] != 0)) # objectif robuste reformule

  ## Objective
  @objective(m, Min, z)

  # Fonction exécutée à chaque callback
  #
  # Input
  # - context_id : permet de déterminer pour quelle raison le callback à
  #   été appelé (solution entière, nouvelle relaxation, ...)
  # - cb_data  :   permet  d'obtenir   diverses  autres   informations
  #   (valeur des bornes inférieures et supérieures , meilleure solution
  #   connue, ...)
  function my_callback_function(cb_data::CPLEX.CallbackContext, context_id::Clong)

      # On teste d'abord si le callback est appelé car une solution entière a été obtenue
      # (cette fonction est définie ci-dessous mais son contenu n'est pas très important)
      if isIntegerPoint(cb_data, context_id)

          # Cette ligne doit être appelée avant de pouvoir récupérer la solution entière
          CPLEX.load_callback_variable_primal(cb_data, context_id)

          # On récupère la valeur de x & y
          x_val = callback_value(cb_data, x)
          y_val = callback_value(cb_data, y)
          z_val = callback_value(cb_data, z)

          #### sous-probleme 1 ####
          # Create the model
          m1 = Model(CPLEX.Optimizer)

          ## Variables
          @variable(m1, delta1[1:n,1:n] >= 0)

          ## Constraints
          @constraint(m1,[i in 1:n, j in 1:n], delta1[i,j] <= D[i,j])
          @constraint(m1, sum(delta1[i,j] for i in 1:n, j in 1:n) <= d1)

          ## Objective
          @objective(m1, Max, sum(d[i,j]*(1 + delta1[i,j])*x_val[i,j] for i in 1:n, j in 1:n if d[i,j] != 0))

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
          @objective(m2, Max, sum((p[i] + delta2[i]*ph[i])*y_val[i] for i in 1:n))

          #resolution
          optimize!(m2)
          z2 = JuMP.objective_value.(m2)
          delta_2 = JuMP.getvalue.( m2[:delta2] )

          if z2 >= S + 1e-6 || z1 > z_val + 1e-6 && z1 < z_val - 1e-6
            con1 = @build_constraint(z >= sum(x_val[i,j]*d[i,j]*(1+ delta_1[i,j]) for i in 1:n, j in 1:n if d[i,j] != 0))
            MOI.submit(m, MOI.LazyConstraint(cb_data), con1)
            con2 = @build_constraint(sum(y_val[v]*(p[v] + ph[v]*delta_2[v]) for v in 1:n) <= S)
            MOI.submit(m, MOI.LazyConstraint(cb_data), con2)
            println("Add constraint x <= 1")
          end
      end
  end

  MOI.set(m, CPLEX.CallbackFunction(), my_callback_function)
  optimize!(m)

end
