#*****************************************************************************
#  Copyright (C) 2017 Lee Worden <worden dot lee at gmail dot com>
#
#  Distributed under the terms of the GNU General Public License (GPL) v.2
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.all import *
import dynamicalsystems
import boxmodel

# function names used by the edge generator: we can't put tuples into
# symbolic expressions directly so we represent them by expressions like
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
def default_compartment_renaming( self, args ):
    #print 'default_compartment_renaming', args
    return tuple(args)

# compartment wrapper converts tuples to a form that can be used in
# symbolic expressions.  By default, we convert the tuple to
# a bm_state function call with the tuple as its arguments.
def default_compartment_wrapper( self, x ):
    return bm_state(*x)

# A mapping from composite compartment names (tuples) into aggregate
# quantities, for the purpose of summing compartments.  This is distinct
# from a compartment renaming function, above.
def default_compartment_aggregation(self, args):
    return tuple(args)

# by default, make tuples of state variables into combined state variables by
# subscripting the first with the others
def default_vertex_namer(self, xx):
    return dynamicalsystems.subscriptedsymbol(*xx)

# sometimes we want to use this instead.  Rather than
# (S,a) -> 'S_a', gives (S,a) -> 'X_{Sa}'.
def X_namer( self, ss ):
    return dynamicalsystems.subscriptedsymbol( 'X', *ss )
x_namer = X_namer

def N_namer( self, ss ):
    return dynamicalsystems.subscriptedsymbol( 'N', *ss )

# We are given an original parameter plus tuples and indices into
# them.  For example, beta, ((S,I), 0), ((I,I), 0), (I,I),
# for the beta parameter for the transition in the first S of (S,I) to
# I as a result of contact with the first I in (I,I).
# We return a bm_param expression that can be used in SR
# (i.e. no tuples in the arguments) and ready to be made into a subscripted
# variable.
# In the simple case, we call this bm_param( beta, I, I ) because
# that's enough to uniquely identify it.  This will become beta_{II}.
# If there are cross interactions
# between for instance the S class from the first model and the I class
# from the second model, we need more detailed subscripting, see below.
def default_param_relabeling( self, p, si, ci, t ):
    #pairwise = iter(vtuple[1:-1])
    #lp = list( zip( pairwise, pairwise ) )
    #print 'simple param relabeling', vtuple, lp
    #pairwise = iter(vtuple[1:-1])
    def allbut( lst, i ):
        return [ x for j,x in enumerate(lst) if j!=i ]
    return bm_param( *(
        [p] +
        ( allbut(si[0],si[1]) if si is not None else [] ) +
        ( allbut(ci[0],ci[1]) if ci is not None else [] )
    ) )
    #return bm_param( *( (vtuple[0],) + reduce( lambda l,m:l+m, (t[:i] + t[i+1:] for t,i in zip(pairwise,pairwise) ) ) ) )

# Fuller parameter subscripting.  For example for
# (beta, (S,I), 1, (I,I), 2, (I,I))
# we return bm_param( beta, S, I, 1, I, I, 2 ).
def full_param_relabeling( self, p, si, ci, t ):
    #print 'full_param_relabeling', vtuple
    from matplotlib.cbook import flatten
    return bm_param( *(flatten( [p, si[0], si[1]] + ([ci[0], ci[1]] if ci is not None else []) )) )

def default_param_namer( self, p_tup ):
    return dynamicalsystems.subscriptedsymbol( *p_tup )

def default_vertex_positioner( self, tuple_graph ):
    """default_vertex_positioner:
    construct an X-Y position for each vertex of the product graph.
    assumes the vertices of the graph are tuples of vertices of the
    models' graphs, in order."""
    directions = [ 0, pi/2, pi/6, pi/3 ] # will add more as needed, I guess
    rotations = [ matrix( [[cos(th),sin(th)],[-sin(th),cos(th)]] ) for th in directions ]
    original_positions = [
	m._graph.get_pos() if m._graph.get_pos() is not None else { v:(i,0) for i,v in enumerate(list(m._sources) + list(m._vars) + list(m._sinks)) }
	for m in self._models
    ]
    #print 'default_vertex_positioner'
    #for m in self._models:
    #    print m._graph.get_pos() if m._graph.get_pos() is not None else 'None'
    #print 'original_positions:',original_positions
    #print 'vertices:', graph.vertices()
    positions = {
	t : sum( r*vector(p[v]) for r,p,v in zip(rotations, original_positions, t.operands()) )
	for t in tuple_graph.vertex_iterator()
    }
    #seq = { v:i for m in product_model._models for i,v in enumerate(m._vars) }
    #positions = {
	#t: [
	#    sum(seq[v]*cos(d) for v,d in zip(t.operands(), directions)),
	#    - sum(seq[v]*sin(d) for v,d in zip(t.operands(), directions))
        #] for t in tuple_graph.vertex_iterator()
    #}
    #print 'positions:', positions
    return positions

def default_sop( self, s_tuple, i, s, t, r ):
    #print 'sop', s_tuple, i, s, t, r
    # return set of t_tuples
    return set( [ s_tuple[:i] + (t,) + s_tuple[i+1:] ] )
def default_bop( self, s_tuple, i, c_tuple, i_, s, t, r ):
    # return set of t_tuples
    return set( [ s_tuple[:i] + (t,) + s_tuple[i+1:] ] )

# inclusive inclusions: interact with the thing whatever position it's in
def tuple_inclusions( self, c, tup, i ):
    #print 'tuple inclusions'
    # slow?
    return [ iota for iota,x in enumerate(tup) if x == c ]

# conservative inclusions: only if it's in the relevant position
def conservative_inclusions( self, c, tup, i ):
    #print 'conservative inclusions'
    return ([i] if tup[i] == c else [])

class BoxModelProductException(Exception): pass

# 'single edge stratifier' is called once for each edge of each
# component model, given a set of product vertices.  It loops
# over those vertices and generates all the product edges involving
# those vertices that are versions of that one component edge.
# yields, for each product edge:
#  source and target vertices, transformed by compartment_renaming,
#   but not wrapped in bm_state or named by vertex_namer.
#  original rate expression together in a tuple with dictionaries
#   specifying vertex replacements, post renaming but pre wrapping and naming,
#   and param replacements.
#   They can't actually be substituted at this stage because they aren't
#   compatible with SR until wrapped.
def default_single_edge_stratifier(
	self, source, target, rate, i, seed_set,
	old_set=set(), cross_interactions=True,
    ):
    # list the compartments and parameters involved in the transition rate
    def getvars(r):
        try: return r.variables()
        except AttributeError: return []
    rate_comps = [ x for x in getvars(rate) if x in self._models[i]._vars ]
    rate_params = set( getvars(rate) ) - set( rate_comps )
    # we can handle constant, linear or bilinear transitions
    if rate_comps == [] or rate_comps == [source]:
	for V in seed_set:
	    s_inclusions = self._inclusions( self, source, V, i )
	    for iota in s_inclusions:
		c_repl = { source: self._compartment_renaming( self, V ) }
	        ws = self._unary_operation( self, V, iota, source, target, rate )
	        for W in ws:
                    Vc = self._compartment_renaming(self, V)
                    Wc = self._compartment_renaming(self, W)
                    if Wc != Vc:
                        # TODO: param_namer
                        p_repl = { p: self._param_relabeling( self, p, (V, iota), None, W ) for p in rate_params }
                        #r = rate.subs( repl )
                        #print r
                        yield ( Vc, Wc, (rate,c_repl,p_repl) )
    elif ( (len(rate_comps) == 2 and source in rate_comps) or
            (len(rate_comps) == 1 and source not in rate_comps) ):
	catalyst, = set(rate_comps) - set([source])
        import itertools
        #print 'edge', rate
	for V,C in itertools.chain( itertools.product(seed_set, old_set), itertools.product(old_set | seed_set, seed_set) ):
            Vc = self._compartment_renaming(self,V)
            Cc = self._compartment_renaming(self,C)
            #print ' source', V, '=>', Vc
            #print ' catalyst', C, '=>', Cc
            ## don't interact with source/sink compartments
            if any( c in m._sources or c in m._sinks for c,m in zip(C,self._models) ): continue
	    # do only the one source inclusion here to avoid duplication
            ## this one line seems to be taking most of the compute time
	    #s_inclusions = set( self._inclusions( self, source, V, i ) ).intersection( set( [i] ) )
            #s_inclusions = ( [i] if i in self._inclusions( self, source, V, i ) else [] )
            s_inclusions = [ z for z in self._inclusions(self, source, V, i) if z==i ]
	    for iota in s_inclusions:
		if cross_interactions:
		    c_inclusions = self._inclusions( self, catalyst, C, i )
		else:
		    c_inclusions = set( self._inclusions( self, catalyst, C, i ) ).intersection( set( [i] ) )
		for iota_ in c_inclusions:
		    c_repl = {
			source: Vc,
			catalyst: Cc
		    }
		    for W in self._binary_operation( self, V, iota, C, iota_, source, target, rate ):
                        Wc = self._compartment_renaming(self, W)
                        #print '  (', W, '=>', Wc, ')'
                        if Wc != Vc:
                            #print '  target', W, '=>', Wc
                            p_repl = { p: self._param_relabeling( self, p, (V, iota), (C, iota_), W ) for p in rate_params }
                            #r = rate.subs( repl )
                            yield ( Vc, Wc, (rate,c_repl,p_repl) )
                            if V == C and iota != iota_:
                                # TODO: is this within-class case right in general?
                                #print '  target again,', Wc
                                p_repl = { p: self._param_relabeling( self, p, (V, iota), (None, iota_), W ) for p in rate_params }
                                yield( Vc, Wc, (rate/catalyst,c_repl,p_repl) )
    else: # wrong variables in rate
	raise BoxModelProductException, "Don't understand rate {0}".format(rate)

def simple_edge_stratifier( self, *args, **kwargs ):
    kwargs['cross_interactions'] = False
    return default_single_edge_stratifier( self, *args, **kwargs )

# This edge generator is called to generate a set of product edges
# given a set of product compartments and the component models.
# It calls its single_edge_generator once for each edge of each
# component model, to generate all the product edges made from
# that original edge.
# returns list of edges in the format yielded by the single edge generator
def default_edge_generator(
	self,
	single_edge_generator=None,
	seed_set=None, cross_interactions=True,
        within_compartment_interactions=False,
    ):
    #print 'default_edge_generator,', ('tuple' if inclusions == tuple_inclusions else 'not tuple'), 'inclusions'
    if single_edge_generator is None:
	single_edge_generator = default_single_edge_stratifier
    # TODO: hacky, fix
    # should be if param_relabeling is none then deduce labeling
    if (
	(self._param_relabeling is default_param_relabeling) and
	cross_interactions and
	any( not set(m1._vars).is_disjoint( set(m2._vars) ) for m1 in self._models for m2 in self._models if m2 is not m1 )
    ):
	self._param_relabeling = full_param_relabeling
    if within_compartment_interactions:
        raise ValueError, 'within_compartment_interactions option not supported in single_edge_generator()'
    import itertools
    def factor_tuples( m ):
        try: return m._raw_tuples
        except AttributeError: return [ (v,) for v in m._vars + list(m._sources)]
    def flatten_tuple( tt ):
        return reduce( lambda x,y:x+y, tt, () )
    if seed_set is None:
	seed_set = set( self._compartment_renaming(self,flatten_tuple(V)) for V in itertools.product( *map(factor_tuples, self._models) ) )
    else:
        seed_set = set( seed_set )
    edges = []
    old_vertices = set()
    print [ (m._sources, m._sinks) for m in self._models ]
    sourcesinks = reduce( lambda x,y:x.union(y), (m._sources | m._sinks for m in self._models), set() )
    while len( seed_set ) > 0:
        #print 'old vertices', old_vertices, 'seed set', seed_set
        # for each edge of each model, we generate a set of derived edges
        # in the product model
        new_edges = list( itertools.chain( *(
	    single_edge_generator(
	        self, v, w, r, i, seed_set,
		old_set=old_vertices,
		cross_interactions=cross_interactions
	    )
	    for i in range(len(self._models))
	    for v,w,r in self._models[i]._graph.edge_iterator()
        ) ) )
	edges += new_edges
	# the edges returned may involve vertices we didn't anticipate
	# so we expand our set of vertices dynamically
	# in which case, we have to do the generation again to include
	# transitions involving the new vertices
        old_vertices |= seed_set
        seed_set = set( v for v,w,r in new_edges ).union( set( w for v,w,r in new_edges ) ) - old_vertices
        seed_set = set( t for t in seed_set if ( set(t) & sourcesinks == set() ) )
	if len(old_vertices) + len(seed_set) > 100 * reduce( lambda a,b:a*b, ( len(m._vars) for m in self._models ) ):
	    raise RuntimeError, 'Recursion produces too many compartments'
    return edges

def simple_edge_generator( self, *args, **kwargs ):
    kwargs['cross_interactions'] = False
    return default_edge_generator( self, *args, **kwargs )

class SubstituteFunctionByName(sage.symbolic.expression_conversions.ExpressionTreeWalker):
    ## dictionary of names of functions, shared across the class
    memoize_fn_names = {}
    def __init__(self, name, replacement):
        """
        Traverse the expression and replace any function that has the given name
        by the given replacement function.
        """
        self._fname = name
        self._replacement_fn = replacement
    def composition(self, ex, operator):
        try: so = self.memoize_fn_names[operator]
        except KeyError, e:
            try: so = operator.name()
            except: so = str(operator)
            self.memoize_fn_names[operator] = so
        #print ' ', so
        if (so == self._fname):
            return self._replacement_fn( *[self(_) for _ in ex.operands()] )
        return super(SubstituteFunctionByName, self).composition(ex,operator)

class CompositeBoxModel(boxmodel.BoxModel):
    """CompositeBoxModel is a boxmodel structure in which compartment names,
    and maybe some parameters, have structure.  Unlike the base BoxModel,
    here we maintain two representations of the graph: one with structure --
    self._tuple_graph, in which each compartment is represented by a 'bm_state'
    function of multiple arguments, and all or some parameters are
    represented by a 'bm_param' function of multiple arguments, and a second
    one, self._graph, the "flow graph", in which those functional forms are collapsed to
    regular variable names.  Calls to bm_state and bm_param are transformed
    to simple SR variable names by the vertex_namer and param_namer functions.
    """
    def __init__(
	    self,
	    tuple_graph,
            raw_tuples,
	    var_tuples = None,
            source_tuples=[],
            sink_tuples=[],
            flow_graph=None,
            vars=None,
            sources=None,
            sinks=None,
            aggregate_names=(),
	    vertex_namer=x_namer,
	    param_namer=default_param_namer,
	    bindings=dynamicalsystems.Bindings()
	):

	self._tuple_graph = tuple_graph
	self._tuple_graph.set_latex_options( edge_labels=True, vertex_shape='rectangle' )

	self._vertex_namer = vertex_namer
	self._param_namer = param_namer

	# for consumption by humans and math software, we map vertex
	# tuples to regular variable names with subscripts
        v_subs = SubstituteFunctionByName( 'bm_state', lambda *x: self._vertex_namer(self, x) )
        p_subs = SubstituteFunctionByName( 'bm_param', lambda *x: self._param_namer(self, x) )
        ## because of suppression of source-sink edges there may be
        ## variables in these three collections that aren't in the
        ## graph

        self._raw_tuples = raw_tuples
        if var_tuples is not None:
	    self._var_tuples = var_tuples
        else:
            self._var_tuples = [ bm_state( *rt ) for rt in self._raw_tuples ]

        vertex_labels = { v: v_subs(v) for v in set( self._var_tuples ).union( source_tuples ).union( sink_tuples ).union( self._tuple_graph.vertex_iterator() ) }
        #print 'vertex labels has', vertex_labels.keys()

        if flow_graph is not None:
            self._graph = flow_graph
        else:
            print 'create substituted graph'
            def bind_edges( edges, bindings ):
                print 'here is bindings:', bindings
                for v, w, e in edges:
                    print e
                    be = bindings(p_subs(v_subs(SR(e))))
                    print be
                    if be != 0: yield ( vertex_labels[v], vertex_labels[w], be )
            flow_edges = list( bind_edges( self._tuple_graph.edge_iterator(), bindings ) )
            #print 'flow edges:', flow_edges
            self._graph = DiGraph(
                flow_edges,
                multiedges=True
            )
            self._graph.set_latex_options( edge_labels=True, vertex_shape='rectangle' )
            positions = self._tuple_graph.get_pos()
            if positions is not None:
                self._graph.set_pos( { vertex_labels[v]:positions[v] for v in vertex_labels.iterkeys() } )

        if vars is not None:
            self._vars = vars
        else:
            self._vars = [ vertex_labels[v] for v in self._var_tuples ]
        self._source_tuples = source_tuples
        if sources is not None:
            self._sources = set( sources )
        else:
            self._sources = set( vertex_labels[v] for v in source_tuples )
        self._sink_tuples = sink_tuples
        if sinks is not None:
            self._sinks = set( sinks )
        else:
            self._sinks = set( vertex_labels[v] for v in sink_tuples )
        print 'vars', self._vars, ' sources', self._sources, 'sinks', self._sinks
        self._aggregate_names = aggregate_names

        print 'extract parameters from rates'
        self.assign_parameters()
	self._bindings = bindings
    def assign_parameters(self):
	# now infer the composite parameters
	# TODO: do this right
	def getvars(r):
	    try: return r.variables()
	    except AttributeError: return []
	self._parameters = reduce( lambda x,y: set(x).union(y), (set(getvars(r)) for f,t,r in self._graph.edges()), set() ).difference( set( self._vars ) )
	#print 'parameters:', self._parameters
    def bind(self, *args, **vargs):
        b = dynamicalsystems.Bindings( *args, **vargs )
        d = deepcopy( self )
        d._graph = DiGraph(
            [ (v,w,b(e)) for v,w,e in d._graph.edge_iterator() ],
            pos = d._graph.get_pos(),
            multiedges=True
        )
        d._bindings = d._bindings + b
        d.assign_parameters()
        return d
    def combine_arrows( self, use_names=() ):
        if use_names: ## if non empty collection of variable names to introduce
            try:
                ## case: use_names is a dictionary { name:sum of things }
                ## transform it to { one of things:name - rest of things }
                ##  to do the substitution that simplifies sum to name
                combine_names_dict = { k:v for k,v in use_names.iteritems() }
                #print '=', use_names
            except AttributeError:
                try:
                    ## case: use_names is a tuple of variables,
                    ## which are keys in self._inclusion_variables
                    use_names[0]
                except TypeError:
                    ## fall back case: use_names is a single variable
                    ## make it into a tuple and continue
                    use_names = ( use_names, )
                #print 'inclusion_variables is', self._inclusion_variables
                combine_names_dict = { v:sum( self._inclusion_variables[v] ) for v in use_names }
            #print 'combined rates:'
            #for v,w,r in b._graph.edge_iterator():
            #    print ' ', r
        else:
            combine_names_dict = {}
        d = {}
        for v,w,r in self._graph.edge_iterator():
            #print r; sys.stdout.flush()
            d[(v,w)] = d.get( (v,w), 0 ) + r
        #print '==', combine_names_dict
        def combine_names(r):
            for k,v in combine_names_dict.iteritems():
                if r == v:
                    r = k
                elif v.operator() == sage.symbolic.operators.add_vararg:
                    rate_vars = set( r.variables() )
                    terms_in_rate = [ x for x in v.operands() if x in rate_vars ]
                    if len(terms_in_rate) > 1:
                        #print ':', r
                        #print ':', k, v
                        expando = { terms_in_rate[0] : k - (v - terms_in_rate[0]) }
                        #print ':', expando
                        r = r.subs( expando ).expand().simplify()
                        #print '::', r
                #else:
                #    r = r.subs( { k:v } )
            return r
        #print 'make ee'; sys.stdout.flush()
        ee = [ (v,w,combine_names(r)) for (v,w),r in d.iteritems() ]
        #print 'make b'; sys.stdout.flush()
        b = boxmodel.BoxModel(
                DiGraph( ee, pos=self._graph.get_pos() ),
                self._vars,
                sources=self._sources,
                sinks=self._sinks,
                aggregate_names=combine_names_dict.keys()
        )
        #print 'done'; sys.stdout.flush()
        return b
	#return self.aggregate_compartments( lambda x:tuple(x), self._param_namer, self._vertex_namer )
    def separate_arrows( self ):
	plus = SR('x+1').operator()
	def arrow_iterator( e ):
	    e = e.expand()
	    if e.operator() == plus:
		for t in e.operands(): yield t
	    else:
		yield e
	return boxmodel.BoxModel( DiGraph(
		[ (v,w,ee) for v,w,e in self._graph.edge_iterator() for ee in arrow_iterator(e) ],
		pos = self._graph._pos,
		multiedges=True
	    ),
	    self._vars
	)
    def aggregate_compartments( self,
            compartment_aggregation=default_compartment_aggregation,
	    param_relabeling_after=lambda self, *x:dynamicalsystems.subscriptedsymbol(*x),
	    vertex_namer_after=default_vertex_namer ):
        aggregate = {}
        for vt in self._tuple_graph.vertex_iterator():
            aggregate.setdefault( self._compartment_aggregation( self, vt.operands() ), [] ).append( vt.operands() )
        ## aggregate is { new_tuple: [old tuples], ... }
        #print 'aggregate:', aggregate
        flow_sums = {}
        for v in self._tuple_graph.vertex_iterator():
	    av = bm_state( *self._compartment_aggregation( self, v.operands() ) )
	    if av not in flow_sums: flow_sums[av] = {}
	    for _,w,e in self._tuple_graph.outgoing_edge_iterator(v):
		aw = bm_state( *self._compartment_aggregation( self, w.operands() ) )
	        flow_sums[av].setdefault( aw, SR(0) )
		#print 'substitute:', e
		es = e.substitute_function( bm_param, lambda *x: param_relabeling_after(self,x) ).substitute_function( bm_state, lambda *x: self._vertex_namer(self,x) ).substitute_function( bm_param, lambda *x: self._param_namer(self,x) )
		#print 'substitute:', e, '|||', es
	        flow_sums[av][aw] += es
	## flow_sums[av][aw] is sum of all transitions from
	## (aggregated tuple) av to aw
	## transitions are in terms of old vertex names and new param names

        ## now do substitutions to transform the transition sums
        agg_eqns, agg_symbols = [], []
        agg_subs = dynamicalsystems.Bindings()
        for newt,oldts in aggregate.items():
            #print 'will combine', sum( self._vertex_namer(*oldt) for oldt in oldts ), '==', vertex_namer_after(*newt)
	    eqn = self._vertex_namer(self, oldts[0]) == vertex_namer_after(self,newt) - sum( self._vertex_namer(self,oldt) for oldt in oldts[1:] )
	    if eqn.lhs() != eqn.rhs():
                agg_symbols.append( self._vertex_namer(self,oldts[0]) )
                agg_eqns.append( eqn )
        agg_graph_dict = {}
        for av, ve in flow_sums.iteritems():
	    vn = av.substitute_function( bm_state, lambda *x:vertex_namer_after(self,x) )
            agg_graph_dict[vn] = {}
            for aw, e in ve.iteritems():
		wn = aw.substitute_function( bm_state, lambda*x:vertex_namer_after(self,x) )
    	        sym = SR.symbol()
                e = self._bindings( e )
    	        #print e
    	        solns = solve( [ sym == e ] + agg_eqns, sym, *agg_symbols, solution_dict=True )
    	        #print 'solve', [ sym == e ] + agg_eqns, ',', [sym] + agg_symbols, '\n  ', solns
    	        if len(solns) == 1:
    	            #print '  ', maxima(sym), [str(k) == str(sym) for k in solns[0].keys()]
    	            el = [ex for k,ex in solns[0].items() if str(k) == str(sym)]
    	            #print '==>', el[0]
    	            agg_graph_dict[vn][wn] = el[0]
		    # TODO: separate transformed sum into multiple arrows
    	        else:
    	            raise RuntimeError, 'Could not simplify expression ' + str(e) + ':' + str(solns)
	#print 'agg_graph_dict', agg_graph_dict
        #self._vc_eqns = vc_eqns
	## make list of transformed variables
	## they are in those dicts, but we want the order
        agg_vars = []
	def xform_state( *t ):
	    s = compartment_aggregation( self, t )
	    return vertex_namer_after( self, s )
        for v in self._var_tuples:
	    av = v.substitute_function( bm_state, xform_state )
	    if av not in agg_vars: agg_vars.append(av)
	#print 'agg_vars', agg_vars
        ## position the aggregates by matching them to a subset of original
	## compartments
	apos = {}
	for t,p in self._tuple_graph.get_pos().iteritems():
	    at = t.substitute_function( bm_state, xform_state )
	    if at not in apos: apos[at] = p
	#print 'apos', apos
	## The transformation produces a flow graph only, not a structured
	## representation
        return boxmodel.BoxModel( DiGraph( agg_graph_dict, pos=apos ), agg_vars )
    def add_transitions( self, trs ):
	# it's an immutable object, so this operation
	# returns a new BoxModel.  trs is a list of (source,target,rate)
	# tuples suitable for adding to self._tuple_graph (rather than _graph).
	fg = deepcopy(self._graph)
	bm_to_vars = lambda e: e.substitute_function( bm_state, lambda *x:self._vertex_namer(self,x)).substitute_function( bm_param, lambda *x: self._param_namer(self,x) )
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
        print 'BoxModelProduct'
	self._models = models
	self._compartment_renaming =  kwargs.pop( 'compartment_renaming', default_compartment_renaming )
        self._compartment_wrapper =   kwargs.pop( 'compartment_wrapper', default_compartment_wrapper )
	self._vertex_namer =          kwargs.pop( 'vertex_namer', default_vertex_namer )
	self._param_relabeling =      kwargs.pop( 'param_relabeling', default_param_relabeling )
        self._param_namer =           kwargs.pop( 'param_namer', default_param_namer )
	self._vertex_positioner =     kwargs.pop( 'vertex_positioner', default_vertex_positioner )
	self._unary_operation =       kwargs.pop( 'unary_operation', default_sop )
	self._binary_operation =      kwargs.pop( 'binary_operation', default_bop )
	self._inclusions =            kwargs.pop( 'inclusions', conservative_inclusions )
	edge_generator =        kwargs.pop( 'edge_generator', default_edge_generator )
	single_edge_generator = kwargs.pop( 'single_edge_generator', None )
        within_compartment_interactions = kwargs.pop( 'within_compartment_interactions', False )
	seed_set =              kwargs.pop( 'seed_set', None )
	if kwargs: raise TypeError, "Unknown named arguments to BoxModelProduct: %s" % str(kwargs)

        print 'edge_generator'
        ## this produces a list of tuples and pre-substitution rate expressions
	raw_edges = list( edge_generator(
	    self,
	    single_edge_generator=single_edge_generator,
            within_compartment_interactions=within_compartment_interactions,
	    seed_set=seed_set
	) )
        #print 'raw edges:\n', raw_edges
        #print [ p_repl for V,W,(r,c_repl,p_repl) in raw_edges ]
        #print [ [ r, r.subs( { c:compartment_wrapper(*C) for c,C in c_repl.iteritems() } ), r.subs( { c:compartment_wrapper(*C) for c,C in c_repl.iteritems() } ).subs( p_repl ) ] for V,W,(r,c_repl,p_repl) in raw_edges ]
        #print [ compartment_wrapper(*V) for V,W,(r,c_repl,p_repl) in raw_edges ]
        #print [ compartment_wrapper(*W) for V,W,(r,c_repl,p_repl) in raw_edges ]
        edges = [ (
            self._compartment_wrapper(self,V),
            self._compartment_wrapper(self,W),
            r.subs( { c:self._compartment_wrapper(self,C) for c,C in c_repl.iteritems() } ).subs( p_repl )
        ) for V,W,(r,c_repl,p_repl) in raw_edges ]
        raw_tuples = list( set( [ V for V,W,rr in raw_edges ] ).union( set( [ W for V,W,rr in raw_edges ] ) ) )

        print 'separate vars, sources, sinks'
	from collections import OrderedDict
	all_vars_d = OrderedDict( (v,None) for v,w,e in edges )
	all_vars_d.update( (w,None) for v,w,e in edges )
        #print all_vars_d.keys()
        vars_d, sources_s, sinks_s = OrderedDict(), set(), set()
        sources = reduce( lambda x,y: x.union( y ), (m._sources for m in self._models), set() )
        sinks = reduce( lambda x,y: x.union( y ), (m._sinks for m in self._models), set() )
        for t in all_vars_d.keys():
            if any( v in sources for v in t.operands() ):
                sources_s = sources_s.union( set( [ t ] ) )
            elif any( v in sinks for v in t.operands() ):
                sinks_s = sinks_s.union( set( [ t ] ) )
            else:
                vars_d[t] = None
        #print vars_d.keys()
        #print sources_s
        #print sinks_s

        ## suppress edges between source/sink compartments
        edges = [ (v,w,e) for v,w,e in edges if v in vars_d or w in vars_d ]

	#print 'edges for product graph:', edges; sys.stdout.flush()
	tuple_graph = DiGraph( edges, multiedges=True )

        print 'super'
	super(BoxModelProduct,self).__init__(
	    tuple_graph=tuple_graph,
            raw_tuples=raw_tuples,
            var_tuples=vars_d.keys(),
            source_tuples=sources_s,
            sink_tuples=sinks_s,
            vertex_namer = self._vertex_namer,
            param_namer = self._param_namer
	)

        print 'vertex_positioner'
	# graphical positions of graph vertices
        #print 'raw tuples', self._raw_tuples
        #print 'tuple graph vertices', self._tuple_graph.vertices()
	self._tuple_graph.set_pos( self._vertex_positioner( self, tuple_graph ) )
        #print 'tuple graph pos', self._tuple_graph.get_pos()
        #print 'graph vertices', self._graph.vertices()
        vsubs = SubstituteFunctionByName( 'bm_state', lambda *x: self._vertex_namer(self,x) )
        self._graph.set_pos( {
            vsubs(vc):p
            for vc,p in self._tuple_graph.get_pos().iteritems()
        } )
        #print 'graph pos', self._graph.get_pos()

	#print 'edges of flow graph:', self._flow_graph.edges(); sys.stdout.flush()
        print 'inclusion bindings'
	self._inclusion_tuples = {}
	for vbm in self._tuple_graph.vertex_iterator():
	    vs = vbm.operands()
	    for v in vs:
		if not v.is_numeric():
	            self._inclusion_tuples.setdefault(v, []).append( vs )

	#print self._inclusion_tuples
	self._inclusion_variables = {
	    k:[self._vertex_namer(self,v) for v in vl]
	    for k,vl in self._inclusion_tuples.iteritems()
	}
        ## don't put them into the bindings by default because of
        ## case where one level of product graph reuses the
        ## compartment names from component graph
	#for k,vl in self._inclusion_variables.iteritems():
        #    print k, 'is sum of', vl
	#    self._bindings.merge_in_place( { k : sum( vl ) } )

        print 'marginal variables'
        self._variable_marginals = {}
        ## given a bm_state() structure, add all its marginals to the dict
        def insert_marginals( *tup ):
            v = self._vertex_namer( self, tup )
            if v not in self._sources and v not in self._sinks:
                for ss in subsets( tup ):
                    if len(ss) > 0:
                        vss = self._vertex_namer( self, ss )
                        if vss != v:
                            self._variable_marginals.setdefault( vss, [] ).append( v )
            return v
        for v in self._tuple_graph.vertex_iterator():
            v.substitute_function( bm_state, insert_marginals )

        print 'marginal parameters'
        self._parameter_marginals = {}
        ## use all the p_repl dicts in the raw edges to list product
        ## parameters that are made from factor parameters
        for V,W,(r,c_repl,p_repl) in raw_edges:
            for fp,pp in p_repl.iteritems():
                ps = pp.substitute_function( bm_param, lambda *x: self._param_namer(self,x) )
                if ps != fp and (fp not in self._parameter_marginals or ps not in self._parameter_marginals[fp]):
                    self._parameter_marginals.setdefault( fp, [] ).append( ps )
        if False:
            ## given a bm_params() structure, add all its marginals to the dict
            import itertools
            def insert_marginals( *psubs ):
                p, subs = ( psubs[0], psubs[1:] )
                param = self._param_namer( self, psubs )
                for ss in subsets( subs ):
                    ## todo: no good
                    pss = self._param_namer( self, tuple(itertools.chain([p],ss)) )
                    if pss != param:
                        self._parameter_marginals.setdefault( pss, [] ).append( param )
                return param
            for v,w,e in self._tuple_graph.edge_iterator():
                e.substitute_function( bm_param, insert_marginals )
        #print self._tuple_graph.edges() # very slow
        print self._parameter_marginals
    def variable_marginals( self, var ):
        if var in self._vars or var in self._sources or var in self._sinks:
            return [ var ]
        else:# if var in self._variable_marginals:
            return self._variable_marginals[var]
        # if not in marginals, raise a KeyError
    def parameter_marginals( self, param ):
        if param in self._parameters:
            return [ param ]
        else:# if param in self._parameter_marginals:
            return self._parameter_marginals[param]
        # if not in marginals, raise a KeyError

def default_sop_strong( self, s_tuple, iset, eis ):
    # return set of t_tuples
    #print 'sop', s_tuple, eis
    tl = list(s_tuple)
    for (v,w,r),i in eis: tl[i] = w
    return set( [ tuple(tl) ] )
def default_bop_strong( self, s_tuple, iset, c_tuple, i_set, eis ):
    # return set of t_tuples
    tl = list(s_tuple)
    for (v,w,r),i in eis: tl[i] = w
    return set( [ tuple(tl) ] )

# this thing is called once for each set of
# component edges, given a set of product vertices.  It loops
# over those vertices and generates all the product edges involving
# those vertices that are versions of that combination of component
# edges.
# TODO: much duplication with the regular single_edge_generator.
# maybe merge
def default_strong_edge_bundle_generator(
	self, eis, seed_set,
	old_set=set(), cross_interactions=True,
        within_compartment_interactions=True
    ):
    if len(eis) == 0: return
    print 'old set', old_set, '/ seed set', seed_set, '/ edges', eis
    # produce all product edges made from these component edges
    # what is the rate of a transition that combines some set of
    # component transitions?
    # who cares? assume this will only be used to combine
    # instances of the same transition, in a power of a single
    # model.
    if len( set( (r for (w,v,r),i in eis) ) ) != 1:
        #print 'too many rates in', eis
        raise RuntimeError, 'overwhelming rate construction problem involving transitions '+str(eis)
    # else: do the replacement in the one rate
    (source,target,rate),i = eis[0]
    # list the compartments involved in the transition
    rate_comps = [ x for x in rate.variables() if x in self._models[i]._vars ]
    rate_params = set( rate.variables() ) - set( rate_comps )
    # we can handle linear or bilinear transitions
    if rate_comps == [source] or rate_comps == []:
	for V in seed_set:
	    if all( i in self._inclusions( self, v, V, i ) for (v,w,r),i in eis ):
                Vc = self._compartment_renaming(self, V)
		c_repl = { s: Vc for s in rate_comps }
	        for W in self._unary_operation( self, V, [i for e,i in eis], eis ):
                    Wc = self._compartment_renaming(self,W)
                    if Wc != Vc:
                        # TODO: param_namer
                        p_repl = { p: self._param_relabeling( self, p, (V, [i for e,i in eis]), None, W ) for p in rate_params }
                        #print V, eis, '=>', W, repl, r.subs(repl)
                        yield ( Vc, Wc, (rate,c_repl,p_repl) )
    elif len(rate_comps) == 2 and source in rate_comps:
	catalyst, = set(rate_comps) - set([source])
        import itertools
	for V,C in itertools.chain( itertools.product(seed_set, old_set), itertools.product(old_set | seed_set, seed_set) ):
	    print 'consider', V, '+', C, 'at', eis
	    # does V have the relevant compartments?
	    # only consider the one inclusion in V, the one given by eis
	    if not all( i in self._inclusions( self, v, V, i ) for (v,w,r),i in eis ):
		print V, 'does not have the inclusions in', eis
	        continue
            Vc = self._compartment_renaming(self,V)
            Cc = self._compartment_renaming(self,C)
	    iota = [ i for e,i in eis ]
	    # in general case, it's all the ways to include the 
	    # rate's catalyst compartment in C
	    # in simple case, those are the inclusions for C as well
	    if not cross_interactions:
                ## simple case: use each level of C at the level where the
                ## edge is given
	        c_inclusions = set( [ tuple( [ i for e,i in eis ] ) ] ) if all( i in self._inclusions( self, catalyst, C, i ) for e,i in eis ) else set()
	        #print catalyst, 'in C:', list( c_inclusions )
	    else:
                ## crossing case: each edge can go at each level where
                ## the catalyst has an inclusion.
		#c_inclusions = Arrangements( inclusions( catalyst, C, i ), len(eis) )
                ## each element of c_inclusions is a tuple iota_ where
                ## iota_[i] is the level at which the catalyst of the ith
                ## edge is included in C
                import itertools
                c_inclusions = set( itertools.product( *(self._inclusions( self, catalyst, C, i ) for e,i in eis ) ) )
	    print 'inclusions of', catalyst, 'in', C, 'for', eis, ':', list( c_inclusions )
	    for iota_ in c_inclusions:
                print 'do edges for', V, iota, C, iota_
    	        c_repl = {
    		    source: Vc,
    		    catalyst: Cc
    	        }
    	        ts = self._binary_operation( self, V, list( iota ), C, list( iota_ ), eis )
		print 'bop returns', ts
		# TODO: check if in-compartment interaction is right
    	        for W in ts:
                    Wc = self._compartment_renaming(self,W)
                    if Wc != Vc:
                        p_repl = { p: self._param_relabeling( self, p, (V, iota), (C, iota_), W ) for p in rate_params }
                        #print V, iota, C, iota_, ':', self._compartment_renaming( self, W )
                        print latex( self._compartment_wrapper(self,Vc) ), '+', latex( self._compartment_wrapper(self,Cc) ), '===>', latex( self._compartment_wrapper(self,Wc) )
                        yield ( Vc, Wc, (rate,c_repl,p_repl) )
                        if within_compartment_interactions and (V == C) and list(iota) != list(iota_):
                            # TODO: is this within-class case right in general?
                            p_repl = { p: self._param_relabeling( p, (V, iota), (None, iota_), W ) for p in rate_params }
                            #print V, iota, iota_, (iota != iota_), '::', self._compartment_renaming( self, W )
                            print latex( self._compartment_wrapper(self,Vc) ), '===>', latex( self._compartment_wrapper(self,Wc) )
                            yield( Vc, Wc, (rate/catalyst,c_repl,p_repl) )
    else: # wrong variables in rate
	raise BoxModelProductException, "Can't stratify rate {0}".format(rate)

def strong_edge_generator(
        self,
	single_edge_generator=None,
	seed_set=None, cross_interactions=True,
        within_compartment_interactions=True
    ):
    if single_edge_generator is None:
	single_edge_generator = default_strong_edge_bundle_generator
    # TODO: hacky, fix
    if (
	(self._param_relabeling is default_param_relabeling) and
	cross_interactions and
	any( not set(m1._vars).is_disjoint( set(m2._vars) ) for m1 in self._models for m2 in self._models if m2 is not m1 )
    ):
	self._param_relabeling = full_param_relabeling
    import itertools
    if seed_set is None:
        #print 'seed set is none'
	seed_set = set( self._compartment_renaming(self,V) for V in ( itertools.product( *((m._vars + list(m._sources)) for m in self._models) ) ) )
    else:
        #print 'seed set is non-none'
        seed_set = set( seed_set )
    #print 'seed set', seed_set
    edges = []
    old_vertices = set()
    import itertools
    sourcesinks = reduce( lambda x,y:x.union(y), (m._sources | m._sinks for m in self._models), set() )
    while len( seed_set ) > 0 :
        # for each edge of each model, we generate a set of derived edges
        # in the product model
        new_edges = list( itertools.chain( *[
	    # call this thing once for each edge in each model, and
	    # in fact more generally once for each tuple ( (e1,i1), (ex,ix), ... )
	    # where ei is an edge in model i for some set of the models
	    single_edge_generator(
	        self, eis, seed_set,
		old_set=old_vertices,
		cross_interactions=cross_interactions,
                within_compartment_interactions=within_compartment_interactions
	    )
	    for iset in map( set, itertools.powerset( set( range(len(self._models)) ) ) )
	    for eis in itertools.product( *([(e,i) for e in self._models[i]._graph.edge_iterator()] for i in iset) )
        ] ) )
	edges += new_edges
        print 'new edges', new_edges
	# the edges returned may involve vertices we didn't anticipate
	# so we expand our set of vertices dynamically
	# in which case, we have to do the generation again to include
	# transitions involving the new vertices
        old_vertices |= seed_set
        seed_set = set( v for v,w,r in new_edges ).union( set( w for v,w,r in new_edges ) ) - old_vertices
        #print 'seed set becomes', seed_set
        # but not vertices made from source and sink compartments
        seed_set = set( [ t for t in seed_set if set(t) & sourcesinks == set() ] )
        #print 'or actually', seed_set
	if len(old_vertices) + len(seed_set) > 100 * reduce( lambda a,b:a*b, ( len(m._vars) for m in self._models ) ):
            #print len(old_vertices) + len(seed_set), ', more than ', 100 * product( len(m._vars) for m in self._models ), ' compartments'; sys.stdout.flush()
	    raise RuntimeError, 'Recursion produces too many compartments'
    return edges
    #return [ ( bm_state( *compartment_renaming( *V ) ), bm_state( *compartment_renaming( *W ) ), r ) for V,W,r in edges ]

def cross( *models ):
    return BoxModelProduct( *models )

def union_edge_generator( self, **whatever ):
    import itertools
    return itertools.chain( *(m._tuple_graph.edge_iterator() for m in self._models) )

def union_positioner( self, graph ):
    return {
	v:(x,-i)
	for i,pos in enumerate( m._graph.get_pos() for m in self._models )
	for v,(x,y) in pos.iteritems()
    }

## bmunion is not a product, but it's easy to implement as one
## it just combines the graphs "side by side"
## (except if they share vertices and possibly edges, they'll be combined)
## note if we didn't combine those it would be direct sum rather than union
def bmunion( *models ):
    """union of the models' vertex and edge sets"""
    return BoxModelProduct( *models,
	vertex_namer = x_namer,
	edge_generator = union_edge_generator,
	vertex_positioner = union_positioner
    )

## TODO: rewrite using new strong_edge_generator with bop
def strong_product( *models, **kwargs ):
    inclusions =        kwargs.pop( 'inclusions', conservative_inclusions )
    within_compartment_interactions = kwargs.pop( 'within_compartment_interactions', True )
    seed_set =          kwargs.pop( 'seed_set', None )
    compartment_renaming =  kwargs.pop( 'compartment_renaming',  default_compartment_renaming )
    compartment_wrapper = kwargs.pop( 'compartment_wrapper', default_compartment_wrapper )
    vertex_namer =      kwargs.pop( 'vertex_namer',      default_vertex_namer )
    param_namer =       kwargs.pop( 'param_namer',       default_param_namer )
    param_relabeling =  kwargs.pop( 'param_relabeling',  full_param_relabeling )
    vertex_positioner = kwargs.pop( 'vertex_positioner', default_vertex_positioner )
    unary_operation =   kwargs.pop( 'unary_operation',   default_sop_strong )
    binary_operation =  kwargs.pop( 'binary_operation',  default_bop_strong )
    single_edge_generator = kwargs.pop( 'single_edge_generator', default_strong_edge_bundle_generator )
    if kwargs: raise TypeError, "Unknown named arguments to strong_product: %s" % str(kwargs)
    return BoxModelProduct( *models,
        inclusions = inclusions,
        within_compartment_interactions = within_compartment_interactions,
        seed_set = seed_set,
	edge_generator = strong_edge_generator,
        single_edge_generator = single_edge_generator,
	compartment_renaming = compartment_renaming,
        compartment_wrapper = compartment_wrapper,
	vertex_namer = vertex_namer,
	param_namer = param_namer,
        param_relabeling = param_relabeling,
	vertex_positioner = vertex_positioner,
        unary_operation = unary_operation,
        binary_operation = binary_operation
    )

def power( model, i, compartment_renaming=lambda self,x:x, param_relabeling=default_param_relabeling, vertex_namer=x_namer ):
    return BoxModelProduct(
	*([model] * i),
	vertex_namer = vertex_namer,
	compartment_renaming=compartment_renaming,
	param_relabeling=param_relabeling
    )

def bmpower( model, i, compartment_renaming=lambda self,x:x, param_relabeling=default_param_relabeling, vertex_namer=x_namer ):
    return power( model, i, compartment_renaming=compartment_renaming, param_relabeling=param_relabeling, vertex_namer=vertex_namer )

def write_product_formula( M1, M2, M12, tfnm, op=r'\times', size1=(3,1), size2=(1,4), size12=(4,4) ):
    # do I need to do this before opening the file, to get preamble right?
    Mtz = M1.plot_boxes( filename=None, figsize=size1, inline=True )
    from dynamicalsystems import latex_output
    ltx = latex_output( tfnm )
    ltx.write(
        '$\\raisebox{-0.5\\height}{\\hbox{', Mtz, '}}\ ', op,
        '\ \\raisebox{-0.5\\height}{\\hbox{',
        M2.transpose_graph().plot_boxes( filename=None, figsize=size2, inline=True ),
        '}}\ =\ \\raisebox{-0.5\\height}{\\hbox{',
        M12.plot_boxes( filename=None, figsize=size12, inline=True ), '}}$'
    )
    ltx.close()

