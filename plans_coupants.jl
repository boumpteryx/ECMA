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

  ## Constraints
  @constraint(m,sum(x[s,j] for j in 1:n if d[s,j] != 0) == 1)
  @constraint(m,sum(x[i,t] for i in 1:n if d[i,t] != 0) == 1)
  @constraint(m,[v in 1:n; v != s && v != t], sum(x[i,v] for i in 1:n if d[i,v] != 0) == sum(x[v,j] for j in 1:n if d[v,j] != 0))
  @constraint(m,[v in 1:n; v != t], y[v] == sum(x[v,j] for j in 1:n if d[v,j] != 0))
  @constraint(m,[i in 1:n; i != t], y[t] == sum(x[i,t] for j in 1:n if d[j,t] != 0))
  @constraint(m, sum(p[i]*y[i] for i in 1:n) <= S)
  @constraint(m, z > sum(d[i][j]*x[i][j] for i in 1:n, j in 1:n if d[i][j] != 0) ) # objectif robuste reformule

  ## Objective
  @Objective(m, Min, z)

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

          # On récupère la valeur de x # TODO
          x_val = callback_value(cb_data, x)

          # TODO
          if x_val > 1 + 1e-6
              con = @build_constraint(x <= 1)
              MOI.submit(m, MOI.LazyConstraint(cb_data), con)
              println("Add constraint x <= 1")
          end
      end
  end

  MOI.set(m, CPLEX.CallbackFunction(), my_callback_function)
  optimize!(m)

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
end
