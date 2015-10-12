from sage.all import *
import boxmodel
from itertools import product
import dynamicalsystems

# in the future: create a StochasticProcess class, put most of this
# there, and have BoxModel produce one?

def km_states( self, N ):
    if sum( self.ode_flow().values() ) != 0:
        # discretize the state space
        ss = [ vector(l) for l in product( *((i/N for i in range(N+1)) for s in self._vars) ) if sum(l) <= 1 ]
    else:
	print 'Reducing dimension to', tuple(self._vars[:-1])
        ss = [ vector(l + (1 - sum(l),)) for l in product( *((i/N for i in range(N+1)) for s in self._vars[:-1]) ) if sum(l) <= 1 ]
    for s in ss: s.set_immutable()
    return ss

boxmodel.BoxModel.stochastic_states = km_states

def km_state_binding_function( self ):
	return lambda s: dynamicalsystems.Bindings( { sn:sv for sn,sv in zip( self._vars, s ) } )

boxmodel.BoxModel.stochastic_state_binding_function = km_state_binding_function

def bm_kolmogorov_eqns( self, N, callback, *xargs ):
    # to convert states' names to integers for indexing
    state_index = { s:i for i,s in enumerate(self._vars) }
    # before beginning, test for conserved total mass
    km_states = self.stochastic_states(N)
    # to bind state variables to their values in a given state s
    bind_state = self.stochastic_state_binding_function()

    # standard basis vectors
    B = list(km_states[0].parent().basis())
    # do the forward or backward construction
    return callback( self, N, km_states, bind_state, state_index, B, *xargs )

def km_var( x, *args ):
    return SR.symbol( x+'_'+ '_'.join( str(a).replace('/','_') for a in args ),
	latex_name = x+'_{%s}' % (''.join(latex(a) for a in args)) )

def p_var( *args ):
    return km_var( 'p', *args )

def km_pos( states, x ):
    try:
        return { km_var(x,*s) : (s[0],-s[1]) for s in states }
    except:
	return { km_var(x,*s) : (s[0],s[0]) for s in states }

def forward_callback( self, N, km_states, bind_state, state_index, B ):
    # construct the forward equations, as a box model, from the state grid
    # and transitions of the original box model.
    km_graph_dict = {
	# for every s in the discretized grid, and every transition e from v to w,
	# a flow of probability mass from s to s - 1/N v + 1/N w
        # at rate e p(s)
	p_var( *s ) : {
	    p_var( *(s - B[ state_index[v] ] / N + B[ state_index[w] ] / N) ) :
		p_var( *s ) * bind_state(s)(e)
		for v,w,e in self._flow_graph.edge_iterator()
		if bind_state(s)(e) != 0
	} for s in km_states
    }
    return boxmodel.BoxModel( DiGraph( km_graph_dict, pos=km_pos( km_states, 'p' ) ),
	[ p_var( *s ) for s in km_states ] )

boxmodel.BoxModel.forward_boxmodel = lambda self, N: bm_kolmogorov_eqns( self, N, forward_callback )

def q_var( *args ):
    return km_var( 'q', *args )

def backward_callback( self, N, km_states, bind_state, state_index, B ):
    # construct the forward equations, from the state grid
    # and transitions of the original box model.
    km_flow = { q_var( *s ) :
	sum( bind_state( s )(e) *
	    ( q_var( *(s - B[state_index[v]]/N + B[state_index[w]]/N) ) - q_var(*s) )
	    for v,w,e in self._flow_graph.edge_iterator()
	    if bind_state(s)(e) != 0 )
	for s in km_states
    }
    return dynamicalsystems.ODESystem(
	km_flow,
	[ q_var( *s ) for s in km_states ]
    )

boxmodel.BoxModel.backward_equations = lambda self, N: bm_kolmogorov_eqns( self, N, backward_callback )
