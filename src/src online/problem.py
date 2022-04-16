#Maintenance planning problem
#https://nbviewer.jupyter.org/github/jckantor/ND-Pyomo-Cookbook/blob/master/notebooks/04.03-Maintenance-Planning.ipynb

import pyomo.environ as pyo
import numpy as np

# problem parameters
T = 10        # planning period from 1..T
M = 2         # length of maintenance period
P = 1         # number of maintenance periods

# daily profits
np.random.seed(1)
c = {k:np.random.uniform() for k in range(1, T+1)}

def maintenance_planning_bigm(c, T, M, P):
    m = pyo.ConcreteModel()

    m.T = pyo.RangeSet(1, T)
    m.Y = pyo.RangeSet(1, T - M + 1)
    m.S = pyo.RangeSet(0, M - 1)

    m.c = pyo.Param(m.T, initialize = c)
    m.x = pyo.Var(m.T, domain=pyo.Binary)
    m.y = pyo.Var(m.T, domain=pyo.Binary)

    # objective
    m.profit = pyo.Objective(expr = sum(m.c[t]*m.x[t] for t in m.T), sense=pyo.maximize)

    # required number P of maintenance starts
    m.sumy = pyo.Constraint(expr = sum(m.y[t] for t in m.Y) == P)

    # no more than one maintenance start in the period of length M
    m.sprd = pyo.Constraint(m.Y, rule = lambda m, t: sum(m.y[t+s] for s in m.S) <= 1)

    # disjunctive constraints
    m.bigm = pyo.Constraint(m.Y, rule = lambda m, t: sum(m.x[t+s] for s in m.S) <= M*(1 - m.y[t]))

    return m

m = maintenance_planning_bigm(c, T, M, P)

#######################################################################################################
#Create solver and solve the model
#solver = pyo.SolverFactory('scip')
solver = pyo.SolverFactory('scip', executable="/homedirs/local/SCIPOptSuite-7.0.1-Linux/bin/scip")
results = solver.solve(m, keepfiles=True, tee=True, load_solutions=False)

results.write()

if results.solver.termination_condition == pyo.TerminationCondition.optimal:
    m.solutions.load_from(results)

m.pprint()
