#*****************************************************************************
#  Copyright (C) 2017 Lee Worden <worden dot lee at gmail dot com>
#
#  Distributed under the terms of the GNU General Public License (GPL) v.2
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.all import *
import boxmodel
import dynamicalsystems

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
	for s,t,w in self._graph.edge_iterator()
    )
    #H = SR(0)
    #for s,t,w in self._graph.edge_iterator():
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
    return dynamicalsystems.HamiltonianODE( H, x_vars, p_vars, bindings=reduce_bindings )

from . import boxkolmogorov
## these are now implemented by JumpProcess
#boxmodel.BoxModel.hamiltonian = lambda self, p_vars, reduce=False : boxkolmogorov.bm_kolmogorov_eqns( self, 1, hamiltonian_callback, p_vars, reduce )

#boxmodel.BoxModel.hamiltonian_system = lambda self, p_vars, reduce=False: boxkolmogorov.bm_kolmogorov_eqns( self, 1, hamiltonian_system_callback, p_vars, reduce )

## add Lagrangian code to JumpProcess sometime if needed

def lagrangian_callback( self, N, km_states, bind_state, state_index, B, p_vars, reduce=True, return_vars=False ):
    # construct lagrangian function
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
    H0 = sum(
	reduce_bindings(r)
	for s,t,r in self._graph.edge_iterator()
    )
    H1 = sum(
	reduce_bindings(r) * (
	    exp( sum(
		p*b for p,b in zip(
		    p_vars,
		    B[ state_index[t] ] - B[ state_index[s] ]
		)
	    ) )
	)
	for s,t,r in self._graph.edge_iterator()
    )
    L1 = sum(
	reduce_bindings(r) * 
	    sum( p*b for p,b in zip(
		    p_vars,
		    B[ state_index[t] ] - B[ state_index[s] ]
	    ) ) *
	    exp( sum( p*b for p,b in zip(
		p_vars,
		B[ state_index[t] ] - B[ state_index[s] ]
	    ) ) )
	for s,t,r in self._graph.edge_iterator()
    )
    L = H0 - H1 + L1
    #for s,t,w in self._graph.edge_iterator():
	#r = B[ state_index[t] ] - B[ state_index[s] ]
	#w = reduce_bindings(w)
	#print 'r,w:', r, w
	#H += w * (exp( sum( p*b for p,b in zip(p_vars,r) ) ) - 1 )
    if return_vars:
	return L, p_vars, reduce_bindings
    else:
	return L

boxmodel.BoxModel.lagrangian = lambda self, p_vars, reduce=False : boxkolmogorov.bm_kolmogorov_eqns( self, 1, lagrangian_callback, p_vars, reduce )
