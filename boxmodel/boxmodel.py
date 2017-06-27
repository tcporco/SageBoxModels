#*****************************************************************************
#  Copyright (C) 2017 Lee Worden <worden dot lee at gmail dot com>
#
#  Distributed under the terms of the GNU General Public License (GPL) v.2
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import graph_latex_patched
from sage.all import *
import dynamicalsystems

from sage.misc.latex import _latex_file_
#from sage.symbolic.relation import solve
from sage.symbolic.function_factory import function

# constant 'enum' values for use with indexing
class deps:
    index, sumover = range(0,2)

def plot_boxmodel_graph( g, filename=None, inline=False, figsize=(6,6), empty_vertices=(), **options ):
    print 'empty vertices:', empty_vertices
    lopts = {
        'graphic_size': figsize,
        'edge_labels': True,
        'edge_thickness' : 0.02,
        #'edge_fills': True,
        #'edge_color': 'white',
        #'edge_thickness': 0.05
        'vertex_shape': 'rectangle',
        'vertices_empty': { x:True for x in empty_vertices },
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
    xp = ''
    if figsize[0] > 6.75 or figsize[1] > 9:
        latex.add_package_to_preamble_if_available('geometry')
        xp = '\\geometry{papersize={' + str(figsize[0] + 10) + 'cm,' + str(figsize[1] + 20) + 'cm}}\n'
    if inline:
        #LT = '\n\\vspace{24pt}\n' + gl + '\n\\vspace{24pt}\n'
        LT = gl
    else:
        LT = _latex_file_( dynamicalsystems.wrap_latex( gl ), title='', extra_preamble=xp )
    if filename is not None:
        #print 'plot to', filename
        LF = open( filename, 'w' )
        LF.write( LT )
        LF.close()
    return LT

## see BoxModel.plot_boxes() method below
## this is a transformation that supports plotting a box model
## graph using per capita flow rates rather than absolute rates
def per_capita_rates(g):
    def to_per_capita(r,s):
        if s in r.variables(): return r/s
        else:
            print 'Warning: rate ', str(r), 'not converted to per capita'
            return r
    return DiGraph(
        [ (v,w,to_per_capita(e,v)) for v,w,e in g.edge_iterator() ],
        multiedges=True,
        pos = g.get_pos()
    )

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
        aggregate_names=(),
        bindings=dynamicalsystems.Bindings()):
	# we are given a directed graph whose vertex labels are state
	# variables, representing fractions of total population,
	# and whose edge labels are flow rates.
	try:
	    graph.edge_iterator()
	except AttributeError:
            try:
                self.__init__( graph._graph, graph._vars, sources=graph._sources, sinks=graph._sinks, aggregate_names=graph._aggregate_names, bindings=graph._bindings )
                return
            except AttributeError:
	        graph = DiGraph(graph)
	self._graph = graph
	self._graph.set_latex_options( edge_labels=True )
	self._sources = Set( sources )
	self._sinks = Set( sinks )
        self._aggregate_names = aggregate_names
	if vars is None:
	    vars = Set( graph.vertices() ) - self._sources - self._sinks
	self._vars = list(vars)
	def getvars(r):
	    try: return r.variables()
	    except AttributeError: return []
	if parameters is None:
	    # avoid namespace confusion with boxmodelproduct.union
            print 'make parameters'; sys.stdout.flush()
	    parameters = list( reduce( lambda x,y: x.union(y), (Set(getvars(r)) for f,t,r in graph.edges()), Set() ) - self._vars - self._aggregate_names )
            print 'made parameters'; sys.stdout.flush()
	self._parameters = parameters
	#print 'parameters:', parameters
        if False:
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
    def bind(self, *args, **vargs):
	bindings = dynamicalsystems.Bindings( *args, **vargs )
	bound_graph = DiGraph( [
	        (bindings(v),bindings(w),bindings(e)) for v,w,e in self._graph.edge_iterator()
	    ],
	    multiedges=True,
	    pos = { bindings(v):p for v,p in self._graph.get_pos().items() } if self._graph.get_pos() is not None else None
	)
	return BoxModel(
	    bound_graph,
	    vars = [ bindings(v) for v in self._vars ],
            sources = Set( bindings(v) for v in self._sources ),
            sinks = Set( bindings(v) for v in self._sinks ),
	    parameters = [ bindings(p) for p in self._parameters ],
	    parameter_dependencies = {
		bindings(p):[(bindings(d),t) for d,t in pd] for p,pd in self._parameter_dependencies.items()
	    },
            aggregate_names = self._aggregate_names,
	    bindings = self._bindings + bindings
	)
    def add_transitions( self, trs ):
	# We take BoxModel to be an immutable object, so this operation
	# returns a new BoxModel.  trs is a list of (source,target,rate)
	# tuples suitable for adding to self._graph
	print 'add_transitions', trs
	print 'parameters before', self._parameters
	nbm = deepcopy(self)
	nbm._graph.add_edges( trs )
	print self._vars
	for f,t,r in trs:
	    try:
		print r
		print r.variables()
		print Set( r.variables() ).difference( Set( self._vars ) )
		nbm._parameters.update( Set( r.variables() ) - self._vars - self._aggregate_names )
	    except AttributeError: pass
	print 'parameters after', nbm._parameters
	return nbm 
    def reorder_latex_variables( self, ex ):
        #return ex
	# Sage likes to write "I S \beta" in unicode or whatever order -
	# we want "\beta S I", and more generally, first parameters and
	# then compartment names, in a sort of order given by the flow
	# of the transitions.  Here we use left-to-right, top-to-bottom
	# order based on the positions given for compartments.
        # this function returns a sort of pseudo-expression that's only
        # suitable for printing, not for doing math with
	try: self._sorter
	except AttributeError:
	    from collections import defaultdict
	    sort_order_map = dict(
                ## parameters first, Greek letters before Roman
                [ (latex(v),(T,T)) for v in self._parameters for T in [-1e+10 if latex(v)[0] == '\\' or latex(v)[0:2] == '{\\' else -0.9e+10] ] +
                ## then compartment names, in order of the graph layout
	        [ (latex(vv),(pp[0],-pp[1])) for vv,pp in self._graph.get_pos().items() ] +
                ## then any aggregate names
                [ (latex(v),(1e+10,1e+10)) for v in self._aggregate_names ]
	    )
	    # this converter is defined later in this file
	    self._sorter = sort_latex_variables(
                ## parameters then compartments
                sort_order_map,
                ## numbers before anything
                order_numbers_as=(-1e+12,-1e+12),
                ## other expressions just after numbers
                order_unknown_as=(-1e+11,-1e+11)
            )
	#print 'use', self._sorter._map, 'on', latex(ex)
	try: return self._sorter( ex )
	except AttributeError: # ex is not an expression
	    return ex
    def plot_boxes( self, filename=None, inline=False, figsize=(6,6), transform_graph=None, **options ):
	g = self._graph
	## apply the user-supplied transform if any
        ## for example, use transform_graph=per_capita_rates to
        ## plot using per capita rather than absolute flow rates
        if transform_graph is not None:
            g = transform_graph(g)
	## tweak the latex representation of the rates
	g = DiGraph(
	    [ self._vars, [ (v,w,self.reorder_latex_variables(e)) for v,w,e in g.edge_iterator() ] ],
            format='vertices_and_edges',
	    multiedges=True,
	    pos = g.get_pos()
	)
        print 'plot_boxes, sources', self._sources, ', sinks', self._sinks
	return plot_boxmodel_graph( g, filename=filename, inline=inline, figsize=figsize, empty_vertices=self._sources | self._sinks, **options )
    def plot( self, *args, **aargs ):
	def lx(s): return '$%s$'%latex(s)
	lfg = DiGraph(
	    [[lx(s) for s in tup] for tup in self._graph.edge_iterator() ],
	    multiedges=True
	)
	vargs = {
	    'edge_labels' : True,
	    'talk' : True
	}
	if 'pos' not in aargs and self._graph.get_pos() is not None:
	    vargs['pos'] = { lx(v) : p for v,p in self._graph.get_pos().items() }
	vargs.update( aargs )
	#print 'plot vargs:', vargs
	return lfg.plot( *args, **vargs )
    def transpose_graph_in_place( self ):
	self._graph.set_pos( { v:(-y,-x) for v,(x,y) in self._graph.get_pos().iteritems() } )
    def transpose_graph( self ):
	nm = deepcopy( self )
        nm.transpose_graph_in_place()
	return nm
    def aggregate_compartments( self, compartment_aggregation ):
        aggregate = {}
        for vt in self._graph.vertex_iterator():
            ## what if vt is simple and doesn't have operands
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
	#return self.aggregate_compartments( lambda x:x )
        d = {}
        for v,w,r in self._graph.edge_iterator():
            d[(v,w)] = d.get( (v,w), 0 ) + r
        ee = [ (v,w,r) for (v,w),r in d.iteritems() ]
        b = BoxModel( DiGraph( ee, pos=self._graph.get_pos() ), self._vars )
        return b
    def separate_arrows( self ):
	plus = SR('x+1').operator()
	def terms_iterator( e ):
	    e = e.expand()
	    if e.operator() == plus:
		for t in e.operands():
		    for tt in terms_iterator(t):
		        yield t
	    else:
		yield e
	return BoxModel( DiGraph(
		[ (v,w,ee) for v,w,e in self._graph.edge_iterator() for ee in terms_iterator(e) ],
		pos = self._graph.get_pos(),
		multiedges=True
	    ),
	    self._vars
	)
    def jump_process(self):
	try:
            self._jump_process
	except AttributeError:
            #print 'making BoxModel JumpProcess'
            nvars = self._sources | self._sinks
            vars = [ v for v in self._vars if v not in nvars ]
            #print 'vars:',vars
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
                [ (to_r(s,t),rate) for s,t,rate in self._graph.edges() ],
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
	for source, target, rate in self._graph.edge_iterator():
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
	#print 'becomes', expr
	return expr

class sort_latex_variables(sage.symbolic.expression_conversions.ExpressionTreeWalker):
    def __init__(self, sort_order_map, order_numbers_as=-oo, order_unknown_as=oo):
        print 'sort_order_map is', sort_order_map
	self._map = sort_order_map
	self._number_order = order_numbers_as
	self._unknown_order = order_unknown_as
	return super(sort_latex_variables,self).__init__(SR(0))
    def arithmetic(self, ex, operator):
        if operator == (2*SR.symbol('x')).operator():
            print 'reorder latex product of', ex.operands()
            ## sort the factors in a multiplication
	    def keyfn(x):
	        try:
		    return self._map[latex(x)]
	        except KeyError:
		    if x.is_numeric(): return self._number_order
		    else: return self._unknown_order
	    ll = sorted( ex.operands(), key=keyfn )
            minusop = (SR.symbol('x')-1).operator() # it's actually +
            ## special case: a factor is -(x-1) :
            ## we will write that as (1-x)
            ## if there's a factor of -1, look for a subtraction
            rev = [ e for e in ll if e.operator() == minusop ] if -1 in ll else []
            if len( rev ) > 0:
                ## there will only be one -1
                ll = [ e for e in ll if e != -1 ]
                rev = rev[:1]
            ## if there are factors of y^-1
            ## we will put those as y in a denominator
            denom = [ d for d in ll if
                d.operator()==(1/SR.symbol('x')).operator()
                and d.operands()[1] == SR(-1)
            ]
            ## function to render each factor in latex
            def to_lx( ex ):
                ## leave out y^-1 factors, they'll come in as y later
                if ex in denom: return ''
                ## subtractions
                if ex.operator() == minusop:
                    ## if reversed, write backwards
                    if ex in rev:
                        return r'\left({}-{}\right)'.format(latex(-ex.operands()[1]),latex(ex.operands()[0]))
                    ## otherwise, write forwards
                    else:
                        return ''.join( (r'\left(',latex(ex),r'\right)') )
                ## write additions
                if ex.operator() == (SR.symbol('x')+1).operator():
                    return r'\left({}\right)'.format('+'.join(ex.operands()))
                ## if it's a compound symbol, put it in parens
                if ex.is_symbol():
                    lx = latex(ex)
                    lxinner = lx
                    while lxinner[0] == '{' and lxinner[-1] == '}':
                        lxinner = lxinner[1:-1]
                    if len(lxinner) > 1 and '_' not in lxinner and '^' not in lxinner and not( lxinner[0] == '\\' and lxinner[1:].isalpha() ):
                        print 'add () to', lxinner
                        return r'\left({}\right)'.format(lxinner)
                    else: return lx
                ## anything else, use default latex rendering
                return latex(ex)
            ## combine the factors in the numerator
            print ll
            lname = ' '.join(to_lx(v) for v in ll)
            ## if any factors in denominator, combine them and make fraction
            if len(denom) > 0:
                print '/', denom
                lden = ' '.join(to_lx(v.operands()[0]) for v in denom)
                lname = r'\frac{'+lname+'}{'+lden+'}'
            #print latex(ex), ' ==> ', lname
	    Msym = SR.symbol( 'M_{}'.format( ZZ.random_element(1e+10) ), latex_name=lname )
            return Msym
        return super(sort_latex_variables,self).arithmetic(ex,operator)
