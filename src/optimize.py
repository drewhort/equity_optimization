'''
This script contains the different optimisation functions
- kolmpollak
- p median
- piecewise linear
- kolm pollak exact linearisation
- kolm pollak minimize number of stores
'''

# import libraries
from pyscipopt import Model, quicksum, multidict, exp
import numpy as np
import itertools
import math
import logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)


def set_kpcoef(origins, destinations, populations, distances, open_current, alpha, kpcoef): # 0 is kpcoef
    model = Model()
    logger.info('set variables')
    # x_d is binary, 1 if destination d is opened, 0 otherwise
    x = {d: model.addVar(vtype="B") for d in destinations}
    # y_o,d is binary, 1 if destination d is assigned to origin o, 0 otherwise
    y = {i: model.addVar(vtype="B") for i in itertools.product(origins, destinations)}

    logger.info('set constraints')
    # constraint: each origin can only be assigned a single destination
    for o in origins:
        model.addCons(quicksum(y[o, d] for d in destinations) == 1)
    
    # constraint: an origin cannot be assigned an unopen destination
    for d in destinations:
        for o in origins:
            model.addCons(y[o, d]-x[d] <= 0)

    # constraint: which destinations are already open
    for d in open_current:
        model.addCons(x[d] == 1)
        
    # Kolm-Pollak Constraint
    model.addCons(quicksum(populations[o]*y[o,d]*exp(alpha*distances[o,d])
                           for d in destinations for o in origins) - kpcoef <= 0)
    
    # NEW objective: minimize the number of destinations
    logger.info('set objective')
    model.setObjective(quicksum(x[d] for d in destinations), 'minimize')
   
    
    # solve the model
    logger.info('optimizing')
    model.optimize()

    # identify which facilities are opened (i.e., their value = 1)
    new_facilities = np.where([int(round(model.getVal(x[d]))) for d in destinations])[0]
    
######## new print statement
    print(new_facilities)
    return(new_facilities) # no longer returning list of facilities????


def kolmpollak(origins, destinations, populations, distances, open_total, open_current, alpha):
    model = Model()
    logger.info('set variables')
    # x_d is binary, 1 if destination d is opened, 0 otherwise
    x = {d: model.addVar(vtype="B") for d in destinations}
    # y_o,d is binary, 1 if destination d is assigned to origin o, 0 otherwise
    y = {i: model.addVar(vtype="B") for i in itertools.product(origins, destinations)}

    logger.info('set constraints')
    # constraint: each origin can only be assigned a single destination
    for o in origins:
        model.addCons(quicksum(y[o, d] for d in destinations) == 1)
    
    # constraint: an origin cannot be assigned an unopen destination
    for d in destinations:
        for o in origins:
            model.addCons(y[o, d]-x[d] <= 0)

    # constraint: the sum of open destinations should equal the number we want to be open
    model.addCons(quicksum(x[d] for d in destinations) == open_total)

    # constraint: which destinations are already open
    for d in open_current:
        model.addCons(x[d] == 1)  
    # formulating the Kolm-Pollak EDE
    z = {}
    # z = {o: model.addVar(vtype="C", name="z(%s)" % (o)) for o in origins}
    w = {o: model.addVar(vtype="C", name="w(%s)" % (o)) for o in origins}
    for o in origins:
        # evaluate the EDE
        z[o] = quicksum(distances[o, d]*y[o, d] for d in destinations)
        exp_power = alpha*z[o]
        model.addCons((w[o] -populations[o]*exp(exp_power)) == 0)

    # objective: minimise the kolmpollak EDE
    logger.info('set objective')
    model.setObjective(quicksum(w[o] for o in origins), 'minimize')

    # solve the model
    logger.info('optimizing')
    model.optimize()

    # identify which facilities are opened (i.e., their value = 1)
    new_facilities = np.where([int(round(model.getVal(x[d]))) for d in destinations])[0]
    
    return(new_facilities)
    
def piecewise_linear(origins, destinations, populations, distances, open_total, open_current, alpha):
    model = Model()
    logger.info('set variables')
    # x_d is binary, 1 if destination d is opened, 0 otherwise
    x = {d: model.addVar(vtype="B") for d in destinations}
    # y_o,d is binary, 1 if destination d is assigned to origin o, 0 otherwise
    y = {i: model.addVar(vtype="B") for i in itertools.product(origins, destinations)}
    

    logger.info('set constraints')
    # constraint: each origin can only be assigned a single destination
    for o in origins:
        model.addCons(quicksum(y[o, d] for d in destinations) == 1)
    
    # constraint: an origin cannot be assigned an unopen destination
    for d in destinations:
        for o in origins:
            model.addCons(y[o, d]-x[d] <= 0)

    # constraint: the sum of open destinations should equal the number we want to be open
    model.addCons(quicksum(x[d] for d in destinations) == open_total)

    # constraint: which destinations are already open
    for d in open_current:
        model.addCons(x[d] == 1)
  
    # formulating the piecewise linear relaxation

    z, b = {}, {}
    w = {o: model.addVar(vtype="C", name="w(%s)" % (o)) for o in origins}
    for o in origins:
        # evaluate the EDE
        z[o] = quicksum(distances[o, d]*y[o, d] for d in destinations)
        #b_o is the max value x_o can take, for linearization
        b[o] = max(distances[o,d] for d in destinations)
        model.addCons(w[o]-populations[o]*(alpha*z[o]+1) >= 0)
       	model.addCons(w[o]-populations[o]*exp(alpha*b[o]/2)*(alpha*z[o]-alpha*b[o]/2+1) >=0)
        model.addCons(w[o]-populations[o]*exp(alpha*2*b[o]/3)*(alpha*z[o]-alpha*2*b[o]/3+1) >=0)
        model.addCons(w[o]-populations[o]*exp(alpha*b[o])*(alpha*z[o]-alpha*b[o]+1) >= 0)

    # objective: minimise the kolmpollak EDE
    logger.info('set objective')
    model.setObjective(quicksum(w[o] for o in origins), 'minimize')

    # solve the model
    logger.info('optimizing')
    model.optimize()

    # identify which facilities are opened (i.e., their value = 1)
    new_facilities = np.where([int(round(model.getVal(x[d]))) for d in destinations])[0]
    
    return(new_facilities)    

def kp_linear_exact(origins, destinations, populations, distances, open_total, open_current, alpha):
    model = Model()
    # model.setPresolve(pyscipopt.scip.PY_SCIP_PARAMSETTING.OFF)
    logger.info('set variables')
    # x_d is binary, 1 if destination d is opened, 0 otherwise
    x = {d: model.addVar(vtype="B") for d in destinations}
    # y_o,d is binary, 1 if destination d is assigned to origin o, 0 otherwise
    y = {i: model.addVar(vtype="B") for i in itertools.product(origins, destinations)}

    logger.info('set constraints')
    # constraint: each origin can only be assigned a single destination
    for o in origins:
        model.addCons(quicksum(y[o, d] for d in destinations) == 1)
    
    # constraint: an origin cannot be assigned an unopen destination
    for d in destinations:
        for o in origins:
            model.addCons(y[o, d]-x[d] <= 0)

    # constraint: the sum of open destinations should equal the number we want to be open
    model.addCons(quicksum(x[d] for d in destinations) == open_total)

    # constraint: which destinations are already open
    for d in open_current:
        model.addCons(x[d] == 1)
  
    #formulating the Kolm-Pollak EDE
    w = {o: model.addVar(vtype="C", name="w(%s)" % (o)) for o in origins}
    for o in origins:
        #evaluate the EDE
        model.addCons((w[o] -quicksum(populations[o]*y[o,d]*exp(alpha*distances[o,d]) for d in destinations)) == 0)

    # objective: minimise the kolmpollak EDE
    # logger.info('set objective')
    # model.setObjective(quicksum(populations[o]*y[o,d]*exp(alpha*distances[o,d]) for d in destinations for o in origins), 'minimize')

    print(alpha)
    model.setObjective(quicksum(w[o] for o in origins), 'minimize')

    # solve the model
    logger.info('optimizing')
    # model.setPresolve(SCIP_PARAMSETTING.OFF)
    # set_optimizer_attribute(model, "presolving/maxrounds", 0)
    model.optimize()

    # identify which facilities are opened (i.e., their value = 1)
    new_facilities = np.where([int(round(model.getVal(x[d]))) for d in destinations])[0]
    print(model. getObjVal ())
    return(new_facilities)

def pmedian(origins, destinations, populations, distances, open_total, open_current):
    # construct model
    model = Model()
    logger.info('set variables')
    # x_d is binary, 1 if destination d is opened, 0 otherwise
    x = {d: model.addVar(vtype="B") for d in destinations}
    # y_o,d is binary, 1 if destination d is assigned to origin o, 0 otherwise
    y = {i: model.addVar(vtype="B") for i in itertools.product(origins, destinations)}
    
    # objective: minimize the population weighted distance
    logger.info('set objectives')
    # tmp = distances.reset_index().merge(populations.reset_index(), how='right', on='origin').set_index(['origin','destination'])
    # tmp['pwd'] = tmp.distance * tmp.population
    model.setObjective(quicksum(
            populations[o]*distances[o,d]*y[o,d] 
                for d in destinations for o in origins), 'minimize')
    # constraint: each origin can only be assigned a single destination
    logger.info('set constraints')
    for o in origins:
        model.addCons(quicksum(y[o,d] for d in destinations) == 1)
    
    # constraint: an origin cannot be assigned an unopen destination
    for d in destinations:
        for o in origins:
            model.addCons(y[o,d]-x[d] <= 0)
    
    # constraint: the sum of open destinations should equal the number we want to be open
    model.addCons(quicksum(x[d] for d in destinations) == open_total)
    
    # constraint: which destinations are already open
    for d in open_current:
        model.addCons(x[d] == 1)
    
    # solve the model
    logger.info('optimizing')
    model.optimize()
    logger.info('optimization complete')
    # identify which facilities are opened (i.e., their value = 1)
    new_facilities = np.where([int(round(model.getVal(x[d]))) for d in destinations])[0]
    return(new_facilities)

