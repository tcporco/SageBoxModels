# requires: boxmodel.py boxkolmogorov.py
# requires:$(SageDynamics)/dynamicalsystems.py $(SageDynamics)/stochasticdynamics.py 
from sage.all import *
import boxmodel, boxkolmogorov

def embedded_markov_callback( self, N, km_states, bind_state, state_index, B, R ):
    # construct matrix for embedded discrete-time markov chain
    # of box model.
    mtx = zero_matrix( R, len(km_states), sparse=True )
    for s in km_states: s.set_immutable()
    #print km_states
    idx = { s:i for i,s in enumerate(km_states) }
    for s in km_states:
	ps = []
        for v,w,e in self._flow_graph.edge_iterator():
	    rt = self._bindings( bind_state(s)(e) )
	    if rt != 0:
	        s1 = s - B[state_index[v]]/N + B[state_index[w]]/N
	        s1.set_immutable()
		ps.append( (s1, rt) )
	if len(ps) == 0:
	    ps = [ (s,1) ]
	pt = sum( rt for s1,rt in ps )
	for s1,rt in ps:
	    mtx[idx[s1],idx[s]] += rt/pt
    return mtx

boxmodel.BoxModel.embedded_discrete_markov_matrix = lambda self, N, R=QQ: boxkolmogorov.bm_kolmogorov_eqns( self, N, embedded_markov_callback, R )

import os, sys
sys.path.append( os.environ['SageDynamics'] )
import dynamicalsystems, stochasticdynamics
class BoxModelDiscreteSimulation( stochasticdynamics.FiniteDimensionalStochasticDynamics ):
    def __init__(self, boxmodel):
	self._model = boxmodel
	self._trans_fast = [ (v,w,fast_callable(r,vars=boxmodel._vars,domain=float)) for v,w,r in boxmodel._flow_graph.edge_iterator() ]
	self._vars_index = { v:i for i,v in enumerate(boxmodel._vars) }
	super( BoxModelDiscreteSimulation, self ).__init__(
	    boxmodel._vars,
	    bindings = boxmodel._bindings # for use by trajectories
	)
    def solve( self, initial_state, start_time=0, end_time=30, grain=1 ):
	# Note total system size N is implied in initial_state
        try: initial_state = [ initial_state(x) for x in self._vars ]
        except TypeError: pass
	self._N = sum( initial_state )
	self._grain = grain
	#print 'sum of', initial_state, 'is', self._N
	return super( BoxModelDiscreteSimulation, self ).solve( initial_state, start_time, end_time )
    def update_time_and_state( self, t, x ):
	#x_bindings = dynamicalsystems.Bindings( zip(
	#    self._model._vars,
	#    (y/self._N for y in x)
	#) )
	#print 'x_bindings:', x_bindings
	# compute rates of all transitions from this state
	transrates = []
	total_r = 0
	for v,w,e in self._trans_fast:
	    er = e( *x )
	    #print 'rate', e, '->', er
	    # should be a number
	    if er > 0:
		total_r += er
		transrates.append( (v,w,total_r) )
	#print 'transrates:', transrates
	# choose one
	r = RR.random_element(0,total_r)
	# note this raises StopIteration if there is none, which we
	# allow to pass because it indicates we've reached an absorbing state
	v,w,tr = next( ( (v,w,tr) for v,w,tr in transrates if tr > r ) )
	#print 'r:',r
	#print 'chose', v, w
	x[self._vars_index[v]] -= self._grain/self._N
	x[self._vars_index[w]] += self._grain/self._N
	#print x
	import numpy.random
	return ( t + numpy.random.exponential( self._grain/total_r ), x )

