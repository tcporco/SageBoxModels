# requires: boxmodel.py boxkolmogorov.py
# requires:$(SageDynamics)/dynamicalsystems.py $(SageDynamics)/stochasticdynamics.py 
from sage.all import *
import boxmodel, dynamicalsystems, boxkolmogorov

class HamiltonianODE(dynamicalsystems.ODESystem):
    def __init__(self, H, x_vars, p_vars, time_variable=SR.symbol('t'), bindings=dynamicalsystems.Bindings()):
	self._H = H
	self._configuration_vars = x_vars
	self._momentum_vars = p_vars
	super(HamiltonianODE,self).__init__( 
	    dict( { x:diff(H,p) for x,p in zip(x_vars,p_vars) },
		**{ p:-diff(H,x) for x,p in zip(x_vars,p_vars) }
	    ),
	    x_vars + p_vars,
	    time_variable=time_variable,
	    bindings=bindings
	)
	print 'time_var is', self._time_variable

def hamiltonian_callback( self, N, km_states, bind_state, state_index, B, p_vars, reduce=True, return_vars=False ):
    # construct hamiltonian function
    # ignore the discretized states
    if reduce and sum( self.ode_flow().values() ) == 0:
	print 'Reducing dimension to', tuple(self._vars[:-1])
	reduce = True
	reduce_bindings = dynamicalsystems.Bindings( { self._vars[-1]: 1 - sum( self._vars[:-1] ) } )
	if len(p_vars) == len(self._vars):
	    reduce_bindings[p_vars[-1]] = 0 
	#print reduce_bindings
	p_vars = p_vars[:len(self._vars)-1]
    else:
	reduce = False
	reduce_bindings = dynamicalsystems.Bindings()
    H = sum(
	reduce_bindings(w) * (
	    exp( sum(
		p*b for p,b in zip(
		    p_vars,
		    B[ state_index[t] ] - B[ state_index[s] ]
		)
	    ) ) - 1
	)
	for s,t,w in self._flow_graph.edge_iterator()
    )
    #H = SR(0)
    #for s,t,w in self._flow_graph.edge_iterator():
	#r = B[ state_index[t] ] - B[ state_index[s] ]
	#w = reduce_bindings(w)
	#print 'r,w:', r, w
	#H += w * (exp( sum( p*b for p,b in zip(p_vars,r) ) ) - 1 )
    if return_vars:
	return H, p_vars, reduce_bindings
    else:
	return H

def hamiltonian_system_callback( self, N, km_states, bind_state, state_index, B, p_vars, reduce=False ):
    H, p_vars, reduce_bindings = hamiltonian_callback( self, N, km_states, bind_state, state_index, B, p_vars, reduce=reduce, return_vars=True )
    x_vars = self._vars[:len(p_vars)]
    return HamiltonianODE( H, x_vars, p_vars, bindings=reduce_bindings )

boxmodel.BoxModel.hamiltonian = lambda self, p_vars, reduce=False : boxkolmogorov.bm_kolmogorov_eqns( self, 1, hamiltonian_callback, p_vars, reduce )

boxmodel.BoxModel.hamiltonian_system = lambda self, p_vars, reduce=False: boxkolmogorov.bm_kolmogorov_eqns( self, 1, hamiltonian_system_callback, p_vars, reduce )
