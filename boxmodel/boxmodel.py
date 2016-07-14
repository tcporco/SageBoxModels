import graph_latex_patched
from sage.all import *
import dynamicalsystems

from sage.misc.latex import _latex_file_
#from sage.symbolic.relation import solve
from sage.symbolic.function_factory import function

# constant 'enum' values for use with indexing
class deps:
    index, sumover = range(0,2)

class BoxModel(SageObject):
    """Parent class for all kinds of box models.
    Note that since this gets its variables from a graph's vertices,
    rather than from indexers, it can't be used in adaptive dynamics.
    Subclasses that generate their boxes, maybe can.
    """
    def __init__(self, graph,
	vars=None,
	parameters=None,
	parameter_dependencies={},
	sources=(),
	sinks=(),
	flow_graph=None,
        bindings=dynamicalsystems.Bindings()):
	# we are given a directed graph whose vertex labels are state
	# variables, representing fractions of total population,
	# and whose edge labels are flow rates.
	try:
	    graph.edge_iterator()
	except AttributeError:
	    graph = DiGraph(graph)
	self._graph = graph
	self._graph.set_latex_options( edge_labels=True )
	self._graph.set_latex_options( vertex_shape='rectangle' )
	self._sources = set( sources )
	self._sinks = set( sinks )
	if vars is None:
	    vars = list( set( graph.vertices() ) - set( sources ) - set( sinks ) )
	self._vars = vars
	def getvars(r):
	    try: return r.variables()
	    except AttributeError: return []
	if parameters is None:
	    # avoid namespace confusion with boxmodelproduct.union
	    parameters = reduce( lambda x,y: x.union(y), (getvars(r) for f,t,r in graph.edges()), set([]) ).difference( self._vars )
	self._parameters = parameters
	#print 'parameters:', parameters
	self._parameter_dependencies = parameter_dependencies
	for p in self._parameters:
	    if p not in self._parameter_dependencies:
		# infer connections between parameters and compartmentalization
		# for now, simple rule:
		# just connect it to the source variable of its arrow
		# TODO: inference including defined quantities like N
		#print 'infer dependencies for parameter', p
		for v,w,e in self._graph.edges():
		    try: vs = getvars(e)
		    except AttributeError: vs = []
		    if p in vs:
			pd = [ v ]
			#print 'found', p, 'in arrow', e
			#print 'infer dependency on', v
			if p in self._parameter_dependencies and self._parameter_dependencies[p] != pd:
			    #print 'but already inferred', self._parameter_dependencies[p]
			    #print 'dependencies of parameter', p, 'are unclear, inferring no dependencies'
			    pd = []
			self._parameter_dependencies[p] = pd
	for p, pd in self._parameter_dependencies.items():
	    try: [ d[0] for d in pd ]
	    except: self._parameter_dependencies[p] = [ (d,deps.index) for d in pd ]
	#print 'parameter dependencies:', self._parameter_dependencies
	self._bindings = bindings
	if self._graph.get_pos() is None:
	    self._graph.set_pos( { v:(i,0) for i,v in enumerate(self._vars) } )
	if flow_graph is None:
	    flow_graph = self._graph
	self._flow_graph = flow_graph
    def bind(self, *args, **vargs):
	bindings = dynamicalsystems.Bindings( *args, **vargs )
	# TODO: this will not bind cross product correctly
	bound_flow_graph = DiGraph( [
	        (bindings(v),bindings(w),bindings(e)) for v,w,e in self._flow_graph.edge_iterator()
	    ],
	    multiedges=True,
	    pos = { bindings(v):p for v,p in self._flow_graph.get_pos().items() } if self._flow_graph.get_pos() is not None else None
	)
	return BoxModel(
	    bound_flow_graph,
	    vars = [ bindings(v) for v in self._vars ],
	    parameters = [ bindings(p) for p in self._parameters ],
	    parameter_dependencies = {
		bindings(p):[(bindings(d),t) for d,t in pd] for p,pd in self._parameter_dependencies.items()
	    },
	    flow_graph = bound_flow_graph,
	    bindings = self._bindings + bindings
	)
    def add_transitions( self, trs ):
	# We take BoxModel to be an immutable object, so this operation
	# returns a new BoxModel.  trs is a list of (source,target,rate)
	# tuples suitable for adding to self._graph (rather than _flow_graph).
	nbm = deepcopy(self)
	nbm._graph.add_edges( trs )
	if self._flow_graph is not self._graph:
	    nbm._flow_graph.add_edges( trs )
	return nbm 
    def tikz_boxes( self, raw=False, inline=False, figsize=(6,6), transform_graph=lambda x:x, **options ):
	if raw:
	    g = self._graph
	else:
	    g = self._flow_graph
	g = transform_graph(g)
	lopts = {
	    'graphic_size': figsize,
	    'edge_labels': True,
	    'edge_thickness' : 0.02,
	    #'edge_fills': True,
	    #'edge_color': 'white',
	    #'edge_thickness': 0.05
	    'vertices_empty': { x:True for x in self._sources | self._sinks },
	    #'vertex_colors': { x:'white' for x in self._sources | self._sinks },
	    #'vertex_label_colors': { x:'white' for x in self._sources | self._sinks }
	}
	graph_latex_patched.setup_latex_preamble()
	gop = graph_latex_patched.GraphLatex(g)
	if inline:
	    lopts['margins'] = (0.5,0.5,0.5,0.5)
	lopts.update( options )
	#g.set_latex_options( **lopts )
	gop.set_options( **lopts )
	gl = gop.latex()
	if inline:
	    #return '\n\\vspace{24pt}\n' + gl + '\n\\vspace{24pt}\n'
	    return gl
	return _latex_file_( dynamicalsystems.wrap_latex( gl ), title='' )
    def plot_boxes( self, filename=None, raw=False, inline=False, **options ):
	# new Tikz/SVG code
	#print 'plot to', filename
        LF = open( filename, 'w' )
	LT = self.tikz_boxes( raw, inline, **options )
	LF.write( LT )
        LF.close()
	return LT
    def plot( self, *args, **aargs ):
	def lx(s): return '$%s$'%latex(s)
	lfg = DiGraph(
	    [[lx(s) for s in tup] for tup in self._flow_graph.edge_iterator() ],
	    multiedges=True
	)
	vargs = {
	    'edge_labels' : True,
	    'talk' : True
	}
	print 'flow_graph pos:', self._flow_graph.get_pos()
	if 'pos' not in aargs and self._flow_graph.get_pos() is not None:
	    vargs['pos'] = { lx(v) : p for v,p in self._flow_graph.get_pos().items() }
	vargs.update( aargs )
	print 'plot vargs:', vargs
	return lfg.plot( *args, **vargs )
    def transpose_graph_in_place( self ):
	self._graph.set_pos( { v:(-y,-x) for v,(x,y) in self._graph.get_pos().iteritems() } )
	if self._flow_graph is not self._graph:
	    self._flow_graph.set_pos( { v:(-y,-x) for v,(x,y) in self._flow_graph.get_pos().iteritems() } )
    def transpose_graph( self ):
	nm = deepcopy( self )
        nm.transpose_graph_in_place()
	return nm
    def aggregate_compartments( self, compartment_aggregation ):
        aggregate = {}
        for vt in self._graph.vertex_iterator():
            aggregate.setdefault( tuple( compartment_aggregation( vt.operands() ) ), [] ).append( vt.operands() )
        ## aggregate is { new vertex: [old vertices], ... }
        print 'aggregate:', aggregate
        flow_sums = {}
        for v in self._graph.vertex_iterator():
	    av = compartment_aggregation( v )
	    if av not in flow_sums: flow_sums[av] = {}
	    for _,w,e in self._graph.outgoing_edge_iterator(v):
		aw = compartment_aggregation( w )
	        flow_sums[av].setdefault( aw, SR(0) )
	        flow_sums[av][aw] += e
	## flow_sums[av][aw] is sum of all transitions from
	## (aggregated vertex) av to aw
	## transitions are in terms of old vertex names

        ## now do substitutions to transform the transition sums
        agg_eqns, agg_symbols = [], []
        agg_subs = dynamicalsystems.Bindings()
        for newt,oldts in aggregate.items():
            print 'will combine', sum( oldts ), '==', newt
            agg_symbols.append( oldts[0] )
            agg_eqns.append( oldts[0] == newt - sum( oldts[1:] ) )
        agg_graph_dict = {}
        for av, ve in flow_sums.iteritems():
            agg_graph_dict[av] = {}
            for aw, e in ve.iteritems():
    	        sym = SR.symbol()
    	        print e,
    	        solns = solve( [ sym == e ] + agg_eqns, sym, *agg_symbols, solution_dict=True )
    	        #print 'solve', [ sym == e ] + agg_eqns, ',', [sym] + agg_symbols, '\n  ', solns
    	        if len(solns) == 1:
    	            #print '  ', maxima(sym), [str(k) == str(sym) for k in solns[0].keys()]
    	            el = [ex for k,ex in solns[0].items() if str(k) == str(sym)]
    	            print '==>', el[0]
    	            agg_graph_dict[av][aw] = el[0]
    	        else:
    	            raise RuntimeError, 'Could not simplify expression ' + str(e) + ':' + str(solns)
	print 'agg_graph_dict', agg_graph_dict
        #self._vc_eqns = vc_eqns
	## make list of transformed variables
	## they are in those dicts, but we want the order
        agg_vars = []
        for v in self._vars:
	    av = compartment_aggregation( v )
	    if av not in agg_vars: agg_vars.append(av)
	print 'agg_vars', agg_vars
        ## position the aggregates by matching them to a subset of original
	## compartments
	apos = {}
	for t,p in self._graph.get_pos().iteritems():
	    at = compartment_aggregation( t )
	    if at not in apos: apos[at] = p
	print 'apos', apos
        return boxmodel.BoxModel( DiGraph( agg_graph_dict, pos=apos ), agg_vars )
    def combine_arrows( self ):
	return self.aggregate_compartments( lambda x:x )
    def separate_arrows( self ):
	plus = SR('x+1').operator()
	def arrow_iterator( e ):
	    e = e.expand()
	    if e.operator() == plus:
		for t in e.operands():
		    for tt in arrow_iterator(t):
		        yield t
	    else:
		yield e
	return BoxModel( DiGraph(
		[ (v,w,ee) for v,w,e in self._graph.edge_iterator() for ee in arrow_iterator(e) ],
		pos = self._graph.get_pos(),
		multiedges=True
	    ),
	    self._vars
	)
    def jump_process(self):
	try:
		self._jump_process
	except AttributeError:
		vars = list( set( self._graph.vertices() ) - set( self._sources ) - set( self._sinks ) )
		var_index = { v:i for i,v in enumerate(vars) }
		for x in self._sources.union( self._sinks ):
			var_index[x] = None
		def to_r( s, t ):
			r = [ 0 for v in vars ]
			if var_index[s] is not None:
				r[var_index[s]] = -1
			if var_index[t] is not None:
				r[var_index[t]] = 1
			return r
		self._jump_process = dynamicalsystems.JumpProcess(
			vars,
			[ (to_r(s,t),rate) for s,t,rate in self._flow_graph.edges() ],
			bindings=self._bindings
		)
	return self._jump_process
    ## for forward_equations see boxkolmogorov.py
    def backward_equations(self, N, q_name='q'):
	return self.jump_process().backward_equations(N,q_name)
    def generator_matrix( self, N, rate_ring=QQ ):
	return self.jump_process().generator_matrix(N, rate_ring)
    def ode_flow(self):
	return self.jump_process().deterministic_flow()
    def ode(self, time_variable=SR.symbol('t'), bindings=dynamicalsystems.Bindings()):
	return self.jump_process().deterministic_ode(time_variable, bindings)
    def micro_transitions( self ):
	# This could produce micro transitions but it isn't right so far
	# TODO: move this to JumpProcess
	# (in addition to making it work)
	ltx = dynamicalsystems.latex_output_base( dynamicalsystems.write_to_string() )
	lines = []
	for source, target, rate in self._flow_graph.edge_iterator():
	    mu = MakeMicro( self, source )
	    ut = mu( rate )
	    print str(ut); sys.stdout.flush()
	    lines += [ r'  & ' + latex(mu.sigma_fn(SR('x'))) + r'\to' + latex(target)
		+ r' \quad\text{ at rate } '
		+ latex( ut )
	    ]
	ltx.write_align( *lines )
	return ltx._output._str

# useful parent class: expression converter that doesn't
# do anything
from sage.symbolic.expression_conversions import SubstituteFunction
class IdentityConverter(SubstituteFunction):
    def __init__(self):
	pass
    def composition(self, ex, operator):
	# override the parent class's function replacing step
	return operator(*map(self, ex.operands()))

class MakeMicro(IdentityConverter):
    _mul = SR('a*b').operator()
    from sage.symbolic.function_factory import function
    delta_fn = function('delta', latex_name=r'\delta')
    sigma_fn = function('sigma', print_latex_func=lambda self, x:r'\sigma_{%s}' % latex(x))
    bm_sum = function( 'sum', print_latex_func=lambda self, x, s, ex:r'\sum_{%s\in %s}%s' %( latex(x), latex(s), latex(ex) ) )
    bm_indicator = function( 'indicator', print_latex_func=lambda self, ev:r'\mathbb{1}\left(%s\right)' % latex(ev) )
    bm_index_param = function( 'bm_index_param' )
    def __init__(self, model, source):
	self._model = model
	self._source = source
	self._working = False
	self._tags = { s : SR.symbol( 'text'+str(s), latex_name=r'\texttt{%s}'%str(s) ) for s in self._model._vars }
    def __call__(self, ex):
	if self._working:
	    return super(MakeMicro,self).__call__(ex)
	self._working = True
	tx = super(MakeMicro,self).__call__( ex / self._source )
	self._working = False
	return (
	   self.bm_indicator( self.sigma_fn( SR.symbol('x') ) == self._tags[self._source] ) *
	   tx.subs( { s : self.bm_sum( SR.symbol('y'), SR.symbol('X'), 1 / SR('N') * self.bm_indicator( self.sigma_fn( SR('y') ) == self._tags[s] ) ) for s in self._model._vars } )
	)
    def arithmetic(self, ex, operator):
	# do special handling to products of things, before processing the
	# things, to catch inner products
	if operator == self._mul:
	    return self.do_inner_product( *ex.operands() )
	else:
	    return reduce( operator, *map(self, ex.operands()) )
    def symbol(self, s):
	return self.do_inner_product( s ) # just in case
    def do_inner_product(self, *args):
	# leave multiplications as is, except in the case of a
	# parameter dependency marked "sumover": convert that from
	# a regular multiplication to an inner product.
	print 'processing product', args
	margs = list(args)
	sumover = []
	dummy_list = ['y', 'z', 'u', 'v', 'w', 's', 't', 'p', 'q', 'r']
	for p,pd in self._model._parameter_dependencies.items():
	    if p in margs:
		print 'found', p, 'in factors:', args
		if all( d in margs + [self._source] for d,x in pd ):
		    print 'found all of its deps', [d for d,x in pd], 'as well'
		    indices_for_p = []
		    p_times = SR(1)
		    for d,ss in pd:
			if ss == deps.sumover:
			    dummy_var = SR.symbol( dummy_list.pop(0) )
			    indices_for_p.append( dummy_var )
			    sumover.append( dummy_var )
			    print 'will sum over', dummy_var, 'in', d; sys.stdout.flush()
			    margs[margs.index(d)] = 1
			    p_times *= self.bm_indicator( self.sigma_fn( dummy_var ) == self._tags[d] )
			    print 'made it through equality'; sys.stdout.flush()
			elif d == self._source:
			    indices_for_p += [SR('x')]
			else:
			    raise ValueError, 'I am confused about dependence on ' + str(d)
		    index_of_p = margs.index(p)
		    margs[index_of_p] = self.bm_index_param( p, *indices_for_p ) * p_times
		    for dv in reversed(sumover):
			margs[index_of_p] = self.bm_sum( dv, SR.symbol('X'), 1 / SR('N') * margs[index_of_p] )
		    margs[index_of_p] = margs[index_of_p].substitute_function(
			self.bm_index_param,
			lambda *args: dynamicalsystems.subscriptedsymbol( *args )
		    )
		    print margs
		else:
		    raise RuntimeError, (
			"Missing parameter dependencies in expression " +
			str( reduce( self._mul, args ) )
		    )
	expr = reduce( self._mul, margs )
	print 'becomes', expr
	return expr
