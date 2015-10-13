from sage.all import *
import dynamicalsystems, boxmodel

# function names used by the edge generator: we can't put tuples into
# symbolic expressions directly so we represent them by instances like
# bm_state( S, a ), bm_param( beta, S, I ), etc.
from sage.symbolic.function_factory import function
bm_state = function( 'bm_state', latex_name='X' )
bm_param = function( 'bm_param', print_latex_func=lambda bmp, p, *xs:'{0}({1})'.format(latex(p),','.join( latex(x) for x in xs )) )

# The general cross product represents its states as tuples of the
# original model's states.  This function can be replaced by another
# one to, for example, make some of these derived compartments synonyms
# for one another, or simply to give them different names.
# This is distinct from a compartment aggregation function, see below,
# which adds compartments together instead of declaring them equivalent.
def default_compartment_renaming(*args):
    return tuple(args)

# A mapping from composite compartment names (tuples) into aggregate
# quantities, for the purpose of summing compartments.  This is distinct
# from a compartment renaming function, above.
def default_compartment_aggregation(*args):
    return tuple(args)

# by default, make tuples of state variables into combined state variables by
# subscripting the first with the others
default_vertex_namer = dynamicalsystems.subscriptedsymbol

# sometimes we want to use this instead.  Rather than
# (S,a) -> 'S_a', gives (S,a) -> 'X_{Sa}'.
def x_namer( *ss ):
    return dynamicalsystems.subscriptedsymbol( 'X', *ss )

# We are given an original parameter plus a list of tuples and indices into
# them.  For example, (beta, (S,I), 0, (I,I), 0, (I,I)),
# for the beta parameter for the transition in the first S of (S,I) to
# I as a result of contact with the first I in (I,I).
# We return a bm_param expression that can be used in SR
# (i.e. no tuples in the arguments) and ready to be made into a subscripted
# variable.
# In the simple case, we call this bm_param( beta, S, I ) because
# that's enough to uniquely identify it.  This will become beta_{SI}.
# If there are cross interactions
# between for instance the S class from the first model and the I class
# from the second model, we need more detailed subscripting, see below.
def default_param_relabeling( *vtuple ):
    print 'simple param relabeling', vtuple
    pairwise = iter(vtuple[1:-1])
    return bm_param( *( (vtuple[0],) + reduce( lambda l,m:l+m, (t[:i] + t[i+1:] for t,i in zip(pairwise,pairwise) ) ) ) )

# Fuller parameter subscripting.  For example for
# (beta, (S,I), 1, (I,I), 2, (I,I))
# we return bm_param( beta, S, I, 1, I, I, 2 ).
def full_param_relabeling( *vtuple ):
    print 'full_param_relabeling', vtuple
    from matplotlib.cbook import flatten
    return bm_param( *(flatten(vtuple[:-1])) )

def default_vertex_positioner( graph, models ):
    """default_vertex_positioner:
    construct an X-Y position for each vertex of the product graph.
    assumes the vertices of the graph are tuples of vertices of the
    models' graphs, in order."""
    directions = [ 0, pi/2, pi/6, pi/3 ] # will add more as needed, I guess
    rotations = [ matrix( [[cos(th),sin(th)],[-sin(th),cos(th)]] ) for th in directions ]
    seq = { v:i for m in models for i,v in enumerate(m._vars) }
    original_positions = [
	m._flow_graph.get_pos() if m._flow_graph.get_pos() is not None else { v:(i,0) for i,v in enumerate(m._vars) }
	for m in models
    ]
    print 'original_positions:',original_positions
    positions = {
	t : sum( r*vector(p[v]) for r,p,v in zip(rotations, original_positions, t.operands()) )
	for t in graph.vertex_iterator()
    }
    #positions = {
	#t: [
	#    sum(seq[v]*cos(d) for v,d in zip(t.operands(), directions)),
	#    - sum(seq[v]*sin(d) for v,d in zip(t.operands(), directions))
        #] for t in graph.vertex_iterator()
    #}
    print 'positions:', positions
    return positions
    import itertools
    positions = {
        tuple(a[1] for a in dtup): [
	    sum(a[0]*cos(d) for a,d in zip(dtup,self.directions)),
	    - sum(a[0]*sin(d) for a,d in zip(dtup,self.directions))
        ] for dtup in itertools.product(
	    *[ enumerate(m._vars) for m in models ]
        )
    }

def default_sop( s_tuple, i, s, t, r ):
    #print 'sop', s_tuple, i, s, t, r
    # return set of t_tuples
    return Set( [ s_tuple[:i] + (t,) + s_tuple[i+1:] ] )
def default_bop( s_tuple, i, c_tuple, i_, s, t, r ):
    # return set of t_tuples
    return Set( [ s_tuple[:i] + (t,) + s_tuple[i+1:] ] )

def tuple_inclusions( c, tup ):
    return [ iota for iota,x in enumerate(tup) if x == c ]

class BoxModelProductException(Exception): pass

# 'single edge stratifier' is called once for each edge of each
# component model, given a set of product vertices.  It loops
# over those vertices and generates all the product edges involving
# those vertices that are versions of that one component edge.
def default_single_edge_stratifier(
	source, target, rate, i, models,
	seed_set, vertex_namer, param_relabeling, compartment_renaming,
	old_set=Set(), cross_interactions=True,
	unary_operation=default_sop, binary_operation=default_bop,
	inclusions=tuple_inclusions
    ):
    print 'old set', old_set, '/ seed set', seed_set
    # produce all product edges made from this component edge
    # list the compartments involved in the transition
    def getvars(r):
        try: return r.variables()
        except AttributeError: return []
    rate_comps = [ x for x in getvars(rate) if x in models[i]._vars ]
    rate_params = set( getvars(rate) ) - set( rate_comps )
    # we can handle constant, linear or bilinear transitions
    if rate_comps == [] or rate_comps == [source]:
	for V in seed_set:
	    s_inclusions = inclusions( source, V )
	    for iota in s_inclusions:
		repl = { source: bm_state( *compartment_renaming( *V ) ) }
	        for W in unary_operation( V, iota, source, target, rate ):
		    print W
		    # TODO: param_namer
		    repl.update( { p: param_relabeling( p, V, iota, W ) for p in rate_params } )
		    yield ( V, W, rate.subs( repl ) )
    elif len(rate_comps) == 2 and source in rate_comps:
	catalyst, = set(rate_comps) - set([source])
        import itertools
	for V,C in itertools.chain( itertools.product(seed_set, old_set), itertools.product(old_set + seed_set, seed_set) ):
	    print 'consider', V, '+', C, 'at', i
	    # do only the one source inclusion here to avoid duplication
	    #s_inclusions = [ iota for iota,x in enumerate(V) if x == source and iota == i ]
	    s_inclusions = Set( inclusions( source, V ) ).intersection( Set( [i] ) )
	    print source, 'in V:', s_inclusions
	    for iota in s_inclusions:
		if cross_interactions:
		    c_inclusions = inclusions( catalyst, C )
		else:
		    c_inclusions = Set( inclusions( catalyst, C ) ).intersection( Set( [i] ) )
		print catalyst, 'in C:', c_inclusions
		for iota_ in c_inclusions:
	            print 'do edges for', V, iota, C, iota_
		    repl = {
			source: bm_state( *compartment_renaming( *V ) ),
			catalyst: bm_state( *compartment_renaming( *C ) )
		    }
		    for W in binary_operation( V, iota, C, iota_, source, target, rate ):
		        repl.update( { p: param_relabeling( p, V, iota, C, iota_, W ) for p in rate_params } )
		        print V, iota, C, iota_, ':', compartment_renaming( *W ), rate.subs( repl )
		        yield ( V, W, rate.subs( repl ) )
		        if V == C:
			    # TODO: is this within-class case right in general?
			    # A: no, needs iota_ somewhere
		            repl.update( { p: param_relabeling( p, V, iota, iota_, W ) for p in rate_params } )
		            print V, iota, iota_, ':', W, rate.subs( repl ) / bm_state(*C)
		            yield( V, W, rate.subs( repl ) / bm_state(*compartment_renaming(*C)) )
    else: # wrong variables in rate
	raise BoxModelProductException, "Don't understand rate {0}".format(rate)

def simple_edge_stratifier( *args, **kwargs ):
    kwargs['cross_interactions'] = False
    return default_single_edge_stratifier( *args, **kwargs )

# This edge generator is called to generate a set of product edges
# given a set of product compartments and the component models.
# It calls its single_edge_generator once for each edge of each
# component model, to generate all the product edges made from
# that original edge.
def default_edge_generator(
	models,
	vertex_namer, param_relabeling, compartment_renaming,
	single_edge_generator=None,
	seed_set=None, cross_interactions=True,
	unary_operation=default_sop, binary_operation=default_bop,
	inclusions=tuple_inclusions
    ):
    if single_edge_generator is None:
	single_edge_generator = default_single_edge_stratifier
    # TODO: hacky, fix
    # should be if param_relabeling is none then deduce labeling
    if (
	(param_relabeling is default_param_relabeling) and
	cross_interactions and
	any( not Set(m1._vars).intersection( Set(m2._vars) ).is_empty() for m1 in models for m2 in models if m2 is not m1 )
    ):
	param_relabeling = full_param_relabeling
    import itertools
    if seed_set is None:
	seed_set = Set( itertools.product( *(m._vars for m in models) ) )
    edges = []
    old_vertices = Set()
    while not seed_set.is_empty():
        # for each edge of each model, we generate a set of derived edges
        # in the product model
        new_edges = list( itertools.chain( *(
	    single_edge_generator(
	        v, w, r, i, models, seed_set,
	        vertex_namer, param_relabeling, compartment_renaming,
		old_set=old_vertices,
	        unary_operation=unary_operation,
		binary_operation=binary_operation,
		cross_interactions=cross_interactions,
		inclusions=inclusions
	    )
	    for i in range(len(models))
	    for v,w,r in models[i]._graph.edge_iterator()
        ) ) )
	edges += new_edges
	# the edges returned may involve vertices we didn't anticipate
	# so we expand our set of vertices dynamically
	# in which case, we have to do the generation again to include
	# transitions involving the new vertices
        old_vertices += seed_set
        seed_set = Set( v for v,w,r in new_edges ).union( Set( w for v,w,r in new_edges ) ) - old_vertices
	if len(old_vertices) + len(seed_set) > 100:
	    raise RuntimeError, 'Recursion produces too many compartments'
    return [ ( bm_state( *compartment_renaming( *V ) ), bm_state( *compartment_renaming( *W ) ), r ) for V,W,r in edges ]

class CompositeBoxModel(boxmodel.BoxModel):
    """CompositeBoxModel is a boxmodel structure in which compartment names,
    and maybe some parameters, have structure.  Unlike the base BoxModel,
    here we maintain two representations of the graph: one with structure --
    in which each compartment is represented by a 'bm_state' function of
    multiple arguments, and all or some parameters are represented by a
    'bm_param' function of multiple arguments.  Those argument lists are
    transformed to simple SR variable names by the vertex_namer and
    param_namer functions."""
    def __init__(
	    self,
	    graph,
	    var_tuples,
	    vertex_namer=x_namer,
	    param_namer=dynamicalsystems.subscriptedsymbol,
	    bindings=dynamicalsystems.Bindings()
	):

	self._graph = graph
	self._graph.set_latex_options( edge_labels=True, vertex_shape='rectangle' )

	# now for consumption by humans and math software, we map vertex
	# tuples to regular variable names with subscripts
	bm_to_vars = lambda e: e.substitute_function( bm_param, param_namer ).substitute_function( bm_state, vertex_namer )
	def bind_edges( edges, bindings ):
	    print 'here is bindings:', bindings
	    for v, w, e in edges:
	        be = bindings(bm_to_vars(e))
		print 'bind', e, 'to', bm_to_vars(e), 'to', be
	        if be != 0: yield ( bm_to_vars(v), bm_to_vars(w), be )
	self._flow_graph = DiGraph(
	    list( bind_edges( self._graph.edge_iterator(), bindings ) ),
	    multiedges=True
	)
	self._flow_graph.set_latex_options( edge_labels=True, vertex_shape='rectangle' )
	positions = self._graph.get_pos()
	if positions is not None:
	    self._flow_graph.set_pos( { bm_to_vars(v):pos for v,pos in positions.items() } )

	self._var_tuples = var_tuples
	self._vars = [ bm_to_vars(v) for v in var_tuples ]
	self._vertex_namer = vertex_namer
	self._param_namer = param_namer
	self._bindings = bindings
    def bind(self, *args, **vargs):
	return CompositeBoxModel(
	    self._graph,
	    self._var_tuples,
	    self._vertex_namer,
	    self._param_namer,
	    self._bindings + dynamicalsystems.Bindings( *args, **vargs )
	)
    def combine_arrows( self ):
	return self.aggregate_compartments( lambda x:tuple(x), self._param_namer, self._vertex_namer )
    def separate_arrows( self ):
	plus = SR(x+1).operator()
	def arrow_iterator( e ):
	    e = e.expand()
	    if e.operator() == plus:
		for t in e.operands(): yield t
	    else:
		yield e
	return BoxModel( DiGraph(
		[ (v,w,ee) for v,w,e in self._flow_graph.edge_iterator() for ee in arrow_iterator(e) ],
		pos = self._flow_graph._pos,
		multiedges=True
	    ),
	    self._vars
	)
    def aggregate_compartments( self, compartment_aggregation=default_compartment_aggregation,
	    param_namer_after=lambda *x:x,
	    vertex_namer_after=default_vertex_namer ):
        aggregate = {}
        for vt in self._graph.vertex_iterator():
            aggregate.setdefault( compartment_aggregation( vt.operands() ), [] ).append( vt.operands() )
        ## aggregate is { new_tuple: [old tuples], ... }
        print 'aggregate:', aggregate
        flow_sums = {}
        for v in self._graph.vertex_iterator():
	    av = bm_state( *compartment_aggregation( v.operands() ) )
	    if av not in flow_sums: flow_sums[av] = {}
	    for _,w,e in self._graph.outgoing_edge_iterator(v):
		aw = bm_state( *compartment_aggregation( w.operands() ) )
	        flow_sums[av].setdefault( aw, SR(0) )
		es = e.substitute_function( bm_param, param_namer_after ).substitute_function( bm_state, self._vertex_namer )
		print 'substitute:', e, '|||', es
	        flow_sums[av][aw] += es
	## flow_sums[av][aw] is sum of all transitions from
	## (aggregated tuple) av to aw
	## transitions are in terms of old vertex names and new param names

        ## now do substitutions to transform the transition sums
        agg_eqns, agg_symbols = [], []
        agg_subs = dynamicalsystems.Bindings()
        for newt,oldts in aggregate.items():
            #print 'will combine', sum( self._vertex_namer(*oldt) for oldt in oldts ), '==', vertex_namer_after(*newt)
	    eqn = self._vertex_namer(*oldts[0]) == vertex_namer_after(*newt) - sum( self._vertex_namer(*oldt) for oldt in oldts[1:] )
	    if eqn.lhs() != eqn.rhs():
                agg_symbols.append( self._vertex_namer(*oldts[0]) )
                agg_eqns.append( eqn )
        agg_graph_dict = {}
        for av, ve in flow_sums.iteritems():
	    vn = av.substitute_function( bm_state, vertex_namer_after )
            agg_graph_dict[vn] = {}
            for aw, e in ve.iteritems():
		wn = aw.substitute_function( bm_state, vertex_namer_after )
    	        sym = SR.symbol()
    	        print e
    	        solns = solve( [ sym == e ] + agg_eqns, sym, *agg_symbols, solution_dict=True )
    	        print 'solve', [ sym == e ] + agg_eqns, ',', [sym] + agg_symbols, '\n  ', solns
    	        if len(solns) == 1:
    	            #print '  ', maxima(sym), [str(k) == str(sym) for k in solns[0].keys()]
    	            el = [ex for k,ex in solns[0].items() if str(k) == str(sym)]
    	            print '==>', el[0]
    	            agg_graph_dict[vn][wn] = el[0]
		    # TODO: separate transformed sum into multiple arrows
    	        else:
    	            raise RuntimeError, 'Could not simplify expression ' + str(e) + ':' + str(solns)
	print 'agg_graph_dict', agg_graph_dict
        #self._vc_eqns = vc_eqns
	## make list of transformed variables
	## they are in those dicts, but we want the order
        agg_vars = []
	def xform_state( *t ):
	    s = compartment_aggregation( t )
	    return vertex_namer_after( *s )
        for v in self._var_tuples:
	    av = v.substitute_function( bm_state, xform_state )
	    if av not in agg_vars: agg_vars.append(av)
	print 'agg_vars', agg_vars
        ## position the aggregates by matching them to a subset of original
	## compartments
	apos = {}
	for t,p in self._graph.get_pos().iteritems():
	    at = t.substitute_function( bm_state, xform_state )
	    if at not in apos: apos[at] = p
	print 'apos', apos
	## The transformation produces a flow_graph only, not a structured
	## representation
        return boxmodel.BoxModel( DiGraph( agg_graph_dict, pos=apos ), agg_vars )
    def add_transitions( self, trs ):
	# it's an immutable object, so this operation
	# returns a new BoxModel.  trs is a list of (source,target,rate)
	# tuples suitable for adding to self._graph (rather than _flow_graph).
	fg = deepcopy(self._flow_graph)
	bm_to_vars = lambda e: e.substitute_function( bm_state, self._vertex_namer).substitute_function( bm_param, self._param_namer )
	fg.add_edges( [ (bm_to_vars(v), bm_to_vars(w), bm_to_vars(e)) for v,w,e in trs ] )
	return boxmodel.BoxModel(
	    fg,
	    vars = self._vars,
	    bindings = self._bindings
        )

class BoxModelProduct(CompositeBoxModel):
    def __init__(
	    self,
	    *models,
	    **kwargs
	):
	self._models = models
	compartment_renaming =  kwargs.pop( 'compartment_renaming', default_compartment_renaming )
	vertex_namer =          kwargs.pop( 'vertex_namer', default_vertex_namer )
	param_relabeling =      kwargs.pop( 'param_relabeling', default_param_relabeling )
	edge_generator =        kwargs.pop( 'edge_generator', default_edge_generator )
	single_edge_generator = kwargs.pop( 'single_edge_generator', None )
	compartment_aggregation = kwargs.pop( 'compartment_aggregation',
	    lambda x:x )
	vertex_positioner =     kwargs.pop( 'vertex_positioner', default_vertex_positioner )
	unary_operation =       kwargs.pop( 'unary_operation', default_sop )
	binary_operation =      kwargs.pop( 'binary_operation', default_bop )
	inclusions =            kwargs.pop( 'inclusions', tuple_inclusions )
	seed_set =              kwargs.pop( 'seed_set', None )
	if kwargs: raise TypeError, "Unknown named arguments to BoxModelProduct: %s" % str(kwargs)

	edges = list( edge_generator(
	    models,
	    vertex_namer,
	    param_relabeling,
	    compartment_renaming,
	    single_edge_generator=single_edge_generator,
	    unary_operation=unary_operation,
	    binary_operation=binary_operation,
	    inclusions=inclusions,
	    seed_set=seed_set
	) )
	print 'edges for product graph:', edges; sys.stdout.flush()
	graph = DiGraph( edges, multiedges=True )

	# graphical positions of graph vertices
	graph.set_pos( vertex_positioner( graph, models ) )

	from collections import OrderedDict
	vars_d = OrderedDict( (v,None) for v,w,e in edges )
	vars_d.update( (w,None) for v,w,e in edges )
	super(BoxModelProduct,self).__init__(
	    graph, vars_d.keys(), vertex_namer
	)

	self._inclusion_tuples = {}
	for vbm in self._graph.vertex_iterator():
	    vs = vbm.operands()
	    for v in vs:
		if not v.is_numeric():
	            self._inclusion_tuples.setdefault(v, []).append( vs )

	#print self._inclusion_tuples
	self._inclusion_variables = {
	    k:[vertex_namer(*v) for v in vl]
	    for k,vl in self._inclusion_tuples.iteritems()
	}
	for k,vl in self._inclusion_variables.iteritems():
	    self._bindings[k] = sum( vl )

	# are vertices from product of graphs?
	representative_tuple = self._graph.vertex_iterator().next().operands()
	if len( representative_tuple ) == len( models ) and all( v in m._vars for v,m in zip( representative_tuple, models ) ):
	    # make list of variables, in order, by taking product
            import itertools
	    varset = set()
	    self._vars = []
	    self._tuples = []
	    for vs in itertools.product( *(m._vars for m in models) ):
	        t = compartment_renaming( *vs )
	        v = vertex_namer( *t )
	        if v not in varset:
		    varset.add(v)
		    self._vars.append(v)
		    self._tuples.append( bm_state( *t ) )
	    print 'made tuples:', self._tuples
	else:
	    # no - the edge generator gave us some other set of vertices
	    self._tuples = [ bm_state(*t) for t in sorted( t.operands() for t in self._graph.vertices()) ]
	    print 'sorted tuples:', self._tuples
	    self._vars = [ vertex_namer(*t.operands()) for t in self._tuples ]

	# now generate all the crossed parameters
	# TODO: do this right
	def getvars(r):
	    try: return r.variables()
	    except AttributeError: return []
	self._parameters = reduce( lambda x,y: set(x).union(y), (getvars(r) for f,t,r in self._flow_graph.edges()), set() ).difference( self._vars )
	#print 'parameters:', self._parameters

	self._vertex_namer = vertex_namer
	self._param_relabeling = param_relabeling

def default_sop_strong( s_tuple, iset, eis ):
    # return set of t_tuples
    print 'sop', s_tuple, eis
    tl = list(s_tuple)
    for (v,w,r),i in eis: tl[i] = w
    return Set( [ tuple(tl) ] )
def default_bop_strong( s_tuple, iset, c_tuple, i_set, eis ):
    # return set of t_tuples
    tl = list(s_tuple)
    for (v,w,r),i in eis: tl[i] = w
    return Set( [ tuple(tl) ] )

# this thing is called once for each set of
# component edges, given a set of product vertices.  It loops
# over those vertices and generates all the product edges involving
# those vertices that are versions of that combination of component
# edges.
# TODO: much duplication with the regular single_edge_generator.
# maybe merge
def default_strong_edge_bundle_generator(
	eis, models,
	seed_set, vertex_namer, param_relabeling, compartment_renaming,
	old_set=Set(), cross_interactions=True,
	unary_operation=default_sop_strong, binary_operation=default_bop_strong,
	inclusions=tuple_inclusions
    ):
    if len(eis) == 0: return
    print 'old set', old_set, '/ seed set', seed_set
    # produce all product edges made from these component edges
    # what is the rate of a transition that combines some set of
    # component transitions?
    # who cares? we assume this will only be used to combine
    # instances of the same transition, in a power of a single
    # model.
    if len( Set( (r for ((w,v,r),i) in eis) ) ) != 1:
        raise RuntimeError, 'overwhelming rate construction problem involving transitions '+str(eis)
    # else: do the replacement in the one rate
    (source,target,rate),i = eis[0]
    # list the compartments involved in the transition
    rate_comps = [ x for x in rate.variables() if x in models[i]._vars ]
    rate_params = set( rate.variables() ) - set( rate_comps )
    # we can handle linear or bilinear transitions
    if rate_comps == [source]:
	for V in seed_set:
	    if all( i in inclusions( v, V ) for (v,w,r),i in eis ):
		repl = { v: bm_state( *compartment_renaming( *V ) ) }
	        for W in unary_operation( V, [i for e,i in eis], eis ):
		    print W
		    # TODO: param_namer
		    repl.update( { p: param_relabeling( p, V, iota, W ) for p in rate_params } )
		    yield ( V, W, r.subs( repl ) )
    elif len(rate_comps) == 2 and source in rate_comps:
	catalyst, = set(rate_comps) - set([source])
        import itertools
	for V,C in itertools.chain( itertools.product(seed_set, old_set), itertools.product(old_set + seed_set, seed_set) ):
	    #print 'consider', V, '+', C, 'at', eis
	    # does V have the relevant compartments?
	    # only consider the one inclusion in V, the one given by eis
	    if not all( i in inclusions( v, V ) for (v,w,r),i in eis ):
		#print V, 'does not have the inclusions in', eis
	        continue
	    iota = [ i for e,i in eis ]
	    # in general case, it's all the ways to include the 
	    # rate's catalyst compartment in C
	    # in simple case, those are the inclusions for C as well
	    if cross_interactions:
		c_inclusions = Arrangements( inclusions( catalyst, C ), len(eis) )
		#print 'inclusions of', eis, 'in', C, ':', list( c_inclusions )
	    else:
	        c_inclusions = Set( tuple( i for e,i in eis ) ) if all( i in inclusions( catalyst, C ) for e,i in eis ) else Set()
	        #print catalyst, 'in C:', list( c_inclusions )
	    for iota_ in c_inclusions:
                #print 'do edges for', V, iota, C, iota_
    	        repl = {
    		    source: bm_state( *compartment_renaming( *V ) ),
    		    catalyst: bm_state( *compartment_renaming( *C ) )
    	        }
    	        ts = binary_operation( V, list( iota ), C, list( iota_ ), eis )
		print 'bop returns', ts
		# TODO: check if in-compartment interaction is right
    	        for W in ts:
    	            repl.update( { p: param_relabeling( p, V, iota, C, iota_, W ) for p in rate_params } )
    	            #print V, iota, C, iota_, ':', compartment_renaming( *W ), rate.subs( repl )
    	            yield ( V, W, rate.subs( repl ) )
    	            if V == C:
    		        # TODO: is this within-class case right in general?
    	                repl.update( { p: param_relabeling( p, V, iota, iota_, W ) for p in rate_params } )
    	                #print V, iota, iota_, ':', W, rate.subs( repl ) / bm_state(*C)
    	                yield( V, W, rate.subs( repl ) / bm_state(*compartment_renaming(*C)) )
    else: # wrong variables in rate
	raise BoxModelProductException, "Don't understand rate {0}".format(rate)

def strong_edge_generator(
	models,
	vertex_namer, param_relabeling, compartment_renaming,
	single_edge_generator=None,
	seed_set=None, cross_interactions=True,
	unary_operation=default_sop, binary_operation=default_bop,
	inclusions=tuple_inclusions
    ):
    if single_edge_generator is None:
	single_edge_generator = default_strong_edge_bundle_generator
    # TODO: hacky, fix
    if (
	(param_relabeling is default_param_relabeling) and
	cross_interactions and
	any( not Set(m1._vars).intersection( Set(m2._vars) ).is_empty() for m1 in models for m2 in models if m2 is not m1 )
    ):
	param_relabeling = full_param_relabeling
    import itertools
    if seed_set is None:
	seed_set = Set( itertools.product( *(m._vars for m in models) ) )
    edges = []
    old_vertices = Set()
    while not seed_set.is_empty():
        # for each edge of each model, we generate a set of derived edges
        # in the product model
        new_edges = list( itertools.chain( *[
	    # call this thing once for each edge in each model, and
	    # in fact more generally once for each tuple ( (e1,i1), (ex,ix), ... )
	    # where ei is an edge in model i for some set of the models
	    single_edge_generator(
	        eis, models, seed_set,
	        vertex_namer, param_relabeling, compartment_renaming,
		old_set=old_vertices,
	        unary_operation=unary_operation,
		binary_operation=binary_operation,
		cross_interactions=cross_interactions,
		inclusions=inclusions
	    )
	    for iset in Subsets( Set( range(len(models)) ) )
	    for eis in itertools.product( *([(e,i) for e in models[i]._graph.edge_iterator()] for i in iset) )
        ] ) )
	edges += new_edges
	# the edges returned may involve vertices we didn't anticipate
	# so we expand our set of vertices dynamically
	# in which case, we have to do the generation again to include
	# transitions involving the new vertices
        old_vertices += seed_set
        seed_set = Set( v for v,w,r in new_edges ).union( Set( w for v,w,r in new_edges ) ) - old_vertices
	if len(old_vertices) + len(seed_set) > 100:
	    raise RuntimeError, 'Recursion produces too many compartments'
    return [ ( bm_state( *compartment_renaming( *V ) ), bm_state( *compartment_renaming( *W ) ), r ) for V,W,r in edges ]

def cross( *models ):
    return BoxModelProduct( *models )

def union_edge_generator( models, vertex_namer, param_namer, compartment_renaming, **whatever ):
    import itertools
    return itertools.chain( *(m._graph.edge_iterator() for m in models) )

def union_positioner( graph, models, compartment_renaming=None ):
    return {
	v:(x,-i)
	for i,pos in enumerate( m._graph.get_pos() for m in models )
	for v,(x,y) in pos.iteritems()
    }

## union is not a product, but it's easy to implement as one
## it just combines the graphs "side by side"
## (except if they share vertices and possibly edges, they'll be combined)
## note if we didn't combine those it would be direct sum rather than union
def union( *models ):
    """union of the models' vertex and edge sets"""
    return BoxModelProduct( *models,
	vertex_namer = x_namer,
	edge_generator = union_edge_generator,
	vertex_positioner = union_positioner
    )

## TODO: rewrite using new strong_edge_generator with bop
def strong_product( *models, **kwargs ):
    compartment_renaming =  kwargs.pop( 'compartment_renaming',  default_compartment_renaming )
    vertex_namer =    kwargs.pop( 'vertex_namer',    default_vertex_namer )
    param_namer =     kwargs.pop( 'param_namer',     default_vertex_namer )
    vertex_positioner = kwargs.pop( 'vertex_positioner', default_vertex_positioner )
    if kwargs: raise TypeError, "Unknown named arguments to BoxModelProduct: %s" % str(kwargs)
    return BoxModelProduct( *models,
	edge_generator = strong_edge_generator,
	vertex_namer = vertex_namer,
	compartment_renaming = compartment_renaming,
	param_namer = param_namer,
	vertex_positioner = vertex_positioner
    )

def power( model, i, compartment_renaming=lambda *x:x, param_relabeling=default_param_relabeling ):
    return BoxModelProduct(
	*([model] * i),
	vertex_namer = x_namer,
	compartment_renaming=compartment_renaming,
	param_relabeling=param_relabeling
    )

def write_product_formula( M1, M2, M12, tfnm, op=r'\times', size1=(3,1), size2=(1,4), size12=(4,4) ):
    # do I need to do this before opening the file, to get preamble right?
    Mtz = M1.tikz_boxes( figsize=size1, inline=True )
    import latex_output
    ltx = latex_output.latex_output( tfnm )
    ltx.write( '$\\raisebox{-0.5\\height}{\\hbox{', Mtz, '}}', op, '\\raisebox{-0.5\\height}{\\hbox{', M2.transpose_graph().tikz_boxes( figsize=size2, inline=True ), '}} = \\raisebox{-0.5\\height}{\\hbox{', M12.tikz_boxes( figsize=size12, inline=True ), '}}$' )
    ltx.close()

