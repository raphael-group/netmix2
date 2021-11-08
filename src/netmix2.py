import numpy as np
from gurobipy import *



def netmix_edgedense(A, rho, resps, alpha=0.05, edge_dense_quad=False, edge_dense_linear=False, output=False, time_limit=None, method=None, thread_count = None, mipgap =None):
    n = A.shape[0]
    clique_size = int(n*alpha)
    print('n: {}, clique_size: {}'.format(n, clique_size))
    
    # set up gurobi
    vertex_inds = [i for i in range(n)]
    inds = vertex_inds

    edge_inds_x, edge_inds_y = np.nonzero(A)
    edge_inds = [(edge_inds_x[t], edge_inds_y[t]) for t in range(len(edge_inds_x)) if edge_inds_x[t] > edge_inds_y[t]]

    # Create a new model
    m = Model("anomaly")
    if not output:
        m.setParam('OutputFlag', 0)
    if time_limit is not None:
        m.setParam('TimeLimit', time_limit)
    if method is not None:
        m.setParam('Method', method)
    if thread_count is not None:
        m.setParam("Threads", thread_count)
    if mipgap is not None:
        m.setParam("MIPGap", mipgap)
    
    print('here1')

    # Create variables
    x=m.addVars(inds, vtype=GRB.BINARY, name="x")
    print('here2')

    # create objective
    w = {i:resps[i] for i in vertex_inds}
    obj_exp = x.prod(w)

    m.setObjective(obj_exp, GRB.MAXIMIZE)
    print('here3')
    
    # size constraint
    m.addConstr( quicksum([x[i] for i in range(n)]) <= clique_size )
    
    # edge density constraint
    if edge_dense_quad:
        LHS = quicksum( [ A[i,j]*x[i]*x[j] for (i,j) in edge_inds ] ) 
        vecsum = quicksum( [x[i] for i in range(n)] )
        RHS = rho * 0.5*( vecsum*vecsum - vecsum )
        m.addConstr( LHS >= RHS  )
    elif edge_dense_linear:
        print('here4')
        LHS = quicksum( [ A[i,j]*x[i]*x[j] for (i,j) in edge_inds ] ) 
        vecsum = quicksum( [x[i] for i in range(n)] )
        RHS = rho * 0.5 * ( vecsum )
        m.addConstr( LHS >= RHS  )
    
    # Optimize model
    m.optimize()
    try:
        estimated_anomaly = [i for i in range(n) if x[i].X > 0.99]
        return estimated_anomaly
    except:
        return []


def netmix_cut(A, rho, resps, alpha=0.05, output=True, time_limit=None):
    n = A.shape[0]
    clique_size = int(n*alpha)
    print('n: {}, clique_size: {}'.format(n, clique_size))
    
    # set up gurobi
    vertex_inds = [i for i in range(n)]
    inds = vertex_inds

    edge_inds_x, edge_inds_y = np.nonzero(A)
    edge_inds = [(edge_inds_x[t], edge_inds_y[t]) for t in range(len(edge_inds_x)) if edge_inds_x[t] > edge_inds_y[t]]

    # Create a new model
    m = Model("anomaly")
    if not output:
        m.setParam('OutputFlag', 0)
    if time_limit is not None:
        m.setParam('TimeLimit', time_limit)
    print('here1')

    # Create variables
    x=m.addVars(inds, vtype=GRB.BINARY, name="x")
    print('here2')

    # create objective
    w = {i:resps[i] for i in vertex_inds}
    obj_exp = x.prod(w)

    m.setObjective(obj_exp, GRB.MAXIMIZE)
    print('here3')
    
    # size constraint
    m.addConstr( quicksum([x[i] for i in range(n)]) <= clique_size )
    
    # cut constraint
    LHS = quicksum([x[i]*x[i] + x[j]*x[j] - 2*x[i]*x[j] for (i,j) in edge_inds])
    RHS=rho*quicksum([x[i] for i in range(n)])
    m.addConstr( LHS <= RHS )
    
    # Optimize model
    m.optimize()
    try:
        estimated_anomaly = [i for i in range(n) if x[i].X > 0.99]
        return estimated_anomaly
    except:
        return []