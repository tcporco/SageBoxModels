#*****************************************************************************
#  Copyright (C) 2017 Lee Worden <worden dot lee at gmail dot com>
#
#  Distributed under the terms of the GNU General Public License (GPL) v.2
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.all import *
from dynamicalsystems import *
from boxmodel import *
from product import *

## inclusions:
##  given compartment 'c' of the i'th factor model,
##  and a compartment 'tup' of the product model,
##  is tup involved in the transitions of c in the product model, and
##  on which level(s).
## practically, will tup be allowed to take the place
##  of c in transitions.
## free_inclusions lets all free model compartments play all roles,
##  though it doesn't cross levels like the pair model does.
def free_inclusions( self, c, tup, i ):
    #print 'free inclusions', c, 'in', tup, ':', [i]
    return [i]

## smarter inclusions for a limited free product that
## only allows sensible sequences of transitions
## it's not trivial - classify each node of the free product
## to a vertex of the Cartesian product graph, and
## allow the arrows that leave that vertex
def classify_vertex( self, t, edge_names={} ):
    try:
        var_key = tuple(tuple(m._vars) for m in self._models)
        edge_names[var_key]
    except KeyError:
        edge_names[var_key] = {
            next( p for p in r.variables() if p not in m._vars ):(v,w,r)
            for m in self._models
            for v,w,r in m._graph.edge_iterator()
        }
    try:
        c = list( t[0] )
    except TypeError: ## shorthand replacing [x] by x
        c = [ t[0] ]
    for eis,cs in t[1:]:
        for e,i in eis:
            v,w,r = edge_names[var_key][e]
            if v == c[i]:
                c[i] = w
            else:
                print "Error in classify_vertex", tup
    return tuple(c)

def cartesian_inclusions( self, c, tup, i ):
    tc = classify_vertex( self, tup )
    return [i] if tc[i] == c else []

def free_uop( self, s_tuple, iset, eis, transition_filter=lambda self,*x:True, edge_names={} ):
    def edge_name( r, exclude=() ):
        if r in edge_names: return edge_names[r]
        pp = [ p for p in r.variables() if p not in exclude ]
        return pp[0]
    s_eis = tuple( (edge_name(r,(v,w)),i) for (v,w,r),i in eis )
    #print 'uop', s_tuple, '+', (s_eis,())
    n_tuple = s_tuple + ( (s_eis,()), )
    if transition_filter( self, (s_tuple, iset), None, eis, n_tuple ):
	#print ': ', n_tuple
	return [ n_tuple ]
    return []

def free_bop( self, s_tuple, iset, c_tuple, i_set, eis, transition_filter=lambda self,*x:True, edge_names={} ):
    def edge_name( r, exclude=() ):
        if r in edge_names: return edge_names[r]
        pp = [ p for p in r.variables() if p not in exclude ]
        return pp[0]
    s_eis = tuple( (edge_name(r,(v,w)),i) for (v,w,r),i in eis )
    #print 'bop', s_tuple, '+', (s_eis,(c_tuple,)) 
    n_tuple = s_tuple + ( (s_eis,(c_tuple,)), )
    ## allow only one kind of transition per arrow
    if ( len( set( e for e,i_ in s_eis ) ) == 1 and
         transition_filter( self, (s_tuple, iset), (c_tuple, i_set), s_eis, n_tuple )
       ):
	#print ': (', n_tuple, ')'
	return [ n_tuple ]
    return []

fw_memo_names={} 
def fw_memoize( tup, v ):
    ## this is not at all elegant, but at least one outside variant
    ## of free_wrapper exists, and calls this, because free_wrapper
    ## stores info here that other functions use
    ## TO DO: put the data in the model where it belongs
    fw_memo_names[tup] = v
def fw_memo( tup ):
    return fw_memo_names[tup]
symbol_by_latex_memo = {}
def sl_memoize( v ):
    symbol_by_latex_memo[latex(v)] = v
def sl_memo( lx ):
    return symbol_by_latex_memo[lx]

## give vertices names by concatenating start compartment with transitions.
## transition names are subscripted according to which of the factor models
## is in play if there are more than one factor model.
## pairwise transitions are marked by the name of the interacting compartment.
def free_wrapper( self, xx ):
    try: v = fw_memo(xx)
    except KeyError:
	## xx[0] is S/I_0/I_1
	## all the rest of xx is (eis,cs) pairs.
	#print 'free wrapper for', xx
	if len(xx) == 1:
	    v = xx[0]
	    #print 'is', v
	else:
	    lname = ''.join(
	        [ latex( xx[0] ) ] +
	        [ ''.join(
		        [ latex( subscriptedsymbol( eis[0][0], *(i for _,i in eis) )
                                if len(self._models > 1)
                                else eis[0][0]
                        ) ] +
		        [ '('+latex(free_wrapper(*c))+')'
                            if len(c) > 1 else
                            ' '+latex(free_wrapper(*c))
                            for c in cs
                        ]
                    )
                    for eis,cs in xx[1:]
	        ]
	    )
            try: v = sl_memo(lname)
            except KeyError:
                v = SR.symbol()
                v = SR.symbol( str(v), latex_name=lname )
                sl_memoize(v)
	    #print 'is', lname
	fw_memoize( xx, v )
    return v

## simplified vertex wrapper that concatenates start compartment with
## transitions but doesn't subscript the transitions
## and doesn't mark interacting compartments.
## to avoid unfortunate duplicate edges, use this with
## rename_boxes_without_catalysts(), below, as compartment renaming.
def free_wrapper_simple( self, tup ):
    try:
        v = fw_memo(tup)
    except KeyError:
        if len(tup) == 1:
            v = tup[0]
        else:
            lname = ''.join( [
                latex( tup[0] )
            ] + [
                latex( eis[0][0] ) for eis,cs in tup[1:]
            ] )
            try: v = sl_memo(lname)
            except:
                v = SR.symbol()
                v = SR.symbol( str(v), latex_name=lname )
                sl_memoize(v)
        fw_memoize(tup, v)
    return v

def rename_boxes_without_catalysts(self, xx):
    return tuple( [ xx[0] ] + [ (eis,()) for eis,cs in xx[1:] ] )

def free_param_relabel( self, p, si, ci, t ):
    return bm_param( *(
        [p] +
        (si[1] + [self._compartment_wrapper(self,si[0])] if si is not None else []) +
        ([ci[0]] if ci is not None and ci[0] is not None else [])
    ) )
    if len(xx) == 4:
        ## unary interaction with source compartment
        p, V, iota, W = xx
        return bm_param( p, *(iota + [self._compartment_wrapper(self,V)]) )
    elif len(xx) == 5:
        ## binary interaction with interaction within source compartment
        p, V, iota, iota_, W = xx
        return bm_param( p, *(iota + [self._compartment_wrapper(self,V)]) )
    elif len(xx) == 6:
        ## binary interaction with source interacting with catalyst
        p, V, iota, C, iota_, W = xx
        return bm_param( p, *(iota + [self._compartment_wrapper(self,V), self._compartment_wrapper(self,C)]) )

p_names = {}
def free_param_namer( self, xx ):
    try: v = p_names[xx]
    except KeyError:
        #print 'param namer', xx
        lname = latex(xx[0])
        subs = [ latex(zz) for zz in xx[1:] ]
	if len(xx) > 1:
            lname = '_'.join( [
		lname,
		( '{' + ''.join(
                    ls if ( len(ls) == 1 or len(subs) == 1 ) else ls.join('()')
                    for ls in subs
                ) + '}' )
	    ] )
        v = SR.symbol()
        v = SR.symbol( str(v), latex_name=lname )
        p_names[xx] = v
        #print 'param name for', xx, 'is', str(v), ':', latex(v)
    return v

## I_01, I_10 used in reductions only
var('I_01, I_10')

memo_pos_base={}
memo_pos={}
memo_pos_places={}
def free_pos( self, graph, init_pos={} ):
    if len(self._models) == 1:
        ## case of a single factor graph is different, simpler
        return free_pos_1( self, graph, init_pos=init_pos )
    ## how to position these compartments.
    ## initial boxes are positioned somewhat smartly
    ## every transition in level 0 adds 1 to x, and
    ## every transition in level 1 adds -1 to y.
    ## overlapping locations are shifted by (0.1,0.1) until unique.
    x_init_pos = {}
    angle = { 0:(0,-1), 1:(1,0), 2:(sqrt(3)/2,1/2), 3:(0,-3) }
    for i,m in enumerate(self._models):
        for v,mp in m._graph.get_pos().iteritems():
            x_init_pos[ (v,) ] = tuple( [ mx*ax for mx,ax in zip(mp,angle[i]) ] )
    x_init_pos.update( init_pos )
    init_pos = x_init_pos
    #print 'init_pos', init_pos
    def free_pos_compartment( xx ):
	try: p = memo_pos[xx]
	except KeyError:
	    #print 'pos', xx
	    if len(xx) == 1:
		#p = { S:(0,0), I_0:(1,0), I_1:(0,-1), I_01:(1,-1), I_10:(1,-1) }[ xx[0] ]
                p = init_pos[ (xx[0],) ]
		side = p[0] + p[1]
	    else:
		yx, eiscs = xx[:-1], xx[-1]
		free_pos_compartment( yx )
		p = vector( memo_pos_base[yx] )
		side = p[1] + p[0]
		#print ' from', yx, 'at', p
		for _,i in eiscs[0]:
                    step = vector( ( (1,0), (0,-1) )[i] )
		    p = p + step
		    if p[0] + p[1] != 0:
                        side = p[0] + p[1]
		    #print ' shift to', p
                p = tuple(p)
	    #print ' base', p
	    memo_pos_base[xx] = p
	    while p in memo_pos_places:
		if side > 0:
	            p = ( p[0] + 0.1, p[1] + 0.15 )
		else:
	            p = ( p[0] - 0.1, p[1] - 0.15 )
	    #print 'is', p
	    memo_pos[xx] = p
	    memo_pos_places[ p ] = 1
	return p
    return { self._compartment_wrapper(*xx) : free_pos_compartment(xx) for xx in self._raw_tuples }

def free_pos_1( self, graph, init_pos={} ):
    ## the factor model's boxes should be laid out horizontally
    ## in the free model they'll be transposed to vertical, and
    ## history will proceed to the right from each.
    yd = { v:x for v,(x,y) in self._models[0]._graph.get_pos().iteritems() }
    pd = {
        self._compartment_wrapper(self, xx) : (len(xx),yd[xx[0]])
        for xx in reversed(self._raw_tuples)
    }
    #print 'free_pos_1 from', pd
    #print 'raw tuples', self._raw_tuples
    ## that dict will have overlapping boxes, they should be spaced out vertically
    dp = { p:[] for p in pd.itervalues() }
    for c,p in pd.iteritems():
        dp[p].append( c )
    for p,cl in dp.iteritems():
        cn = len(cl)
        for i,c in enumerate(sorted(cl,key=lambda v:latex(v))):
            pd[c] = (p[0],p[1] + i*1.0/cn)
    return pd

def limited_uop( st, transition_filter=lambda *x:True ):
    def llu( self, s_tuple, iset, eis ):
        if s_tuple in st: return free_uop( self, s_tuple, iset, eis, transition_filter=transition_filter )
        else: return []
    return llu

def limited_bop( st, transition_filter=lambda *x:True ):
    def llb( self, s_tuple, iset, c_tuple, i_set, eis ):
        if s_tuple in st and c_tuple in st: return free_bop( self, s_tuple, iset, c_tuple, i_set, eis, transition_filter=transition_filter )
        else: return []
    return llb

## for use as a transition filter
def limit_to_iteration( max_iteration, ellipsis_iteration ):
    iteration_dict = {}
    def lim_it( self, si, ci, eis, t_tuple ):
        ## s_tuple, iset, c_tuple, i_set, s_eis, t_tuple
        s_tuple = si[0]
        try:
            s_it = iteration_dict[s_tuple]
        except KeyError:
            s_it = 0
            iteration_dict[s_tuple] = s_it
        if ci is not None and ci[0] is not None:
            c_tuple = ci[0]
            try:
                c_it = iteration_dict[c_tuple]
            except KeyError:
                c_it = 0
                iteration_dict[c_tuple] = c_it
        else:
            c_it = -1
        t_it = max(s_it,c_it) + 1
        iteration_dict[t_tuple] = t_it
        return ( t_it <= max_iteration )
    def ell_it( self, tup ):
        #print 'seek', tup, 'in', iteration_dict
        return ( iteration_dict[tup] == ellipsis_iteration )
    return lim_it, ell_it

def free_product( *models, **kwargs ):
    inclusions =        kwargs.pop( 'inclusions', free_inclusions )
    transition_filter = kwargs.pop( 'transition_filter', None )
    max_iteration =     kwargs.pop( 'max_iteration', None )
    if transition_filter is None:
        if max_iteration is not None:
            ellipsis_iteration = kwargs.pop( 'ellipsis_iteration', None )
            transition_filter, ellipsis_filter = limit_to_iteration( max_iteration, ellipsis_iteration )
        else:
            transition_filter = lambda *x:True
            ellipsis_filter = lambda *x:False
    unary_operation =   kwargs.pop( 'unary_operation', lambda *x: free_uop(*x, transition_filter=transition_filter) )
    binary_operation =  kwargs.pop( 'binary_operation', lambda *x: free_bop(*x, transition_filter=transition_filter) )
    seed_set =          kwargs.pop( 'seed_set', None )
    compartment_renaming = kwargs.pop( 'compartment_renaming', lambda self,x: x )
    compartment_wrapper = kwargs.pop( 'compartment_wrapper', free_wrapper )
    vertex_namer =      kwargs.pop( 'vertex_namer', lambda self,x:x )
    param_relabeling =  kwargs.pop( 'param_relabeling', free_param_relabel )
    param_namer =       kwargs.pop( 'param_namer', free_param_namer )
    vertex_positioner = kwargs.pop( 'vertex_positioner', free_pos )
    if seed_set is not None:
        ## like, in case there's one model and there's a seed set of
        ## bare vertices instead of singleton tuples, convert them
        def mktuple( s ):
            try:
                s[0]
            except TypeError:
                s = (s,)
            return s
        seed_set = [ mktuple(s) for s in seed_set ]
    PM = strong_product( *models,
        inclusions = inclusions,
        unary_operation = unary_operation,
        binary_operation = binary_operation,
	seed_set = seed_set,
        compartment_renaming = compartment_renaming,
        compartment_wrapper = compartment_wrapper,
        vertex_namer = vertex_namer,
        param_relabeling = param_relabeling,
        param_namer = param_namer,
        vertex_positioner = vertex_positioner,
	single_edge_generator = lambda self, eis, *a, **g: default_strong_edge_bundle_generator( self, eis, *a, **g ) if len( Set( (r for (v,w,r),i in eis) ) ) == 1 else []
    )
    PM._ellipsis_tuples = tuple( [ v for v in PM._raw_tuples if ellipsis_filter( PM, v ) ] )
    PM._ellipsis_vars = tuple( [ vertex_namer( compartment_wrapper( PM, v ) ) for v in PM._ellipsis_tuples ] )
    print 'ellipses at', PM._ellipsis_vars
    return PM

## to do: redo for general case
# vertex coloring in pair model allows four colors
def vertex_color( v, memo={} ):
    if v in memo: return memo[v]
    #print 'vertex color of', v
    l = list(v)
    try:
	## I_01, I_10 are generated by reductions.
        color = { S:[S,S], I_0:[I,S], I_1:[S,I], I_01:[I,I], I_10:[I,I] }[l.pop(0)]
    except KeyError:
	raise ValueError, 'Invalid tuple in vertex_color, ' + str(v)
    #print ' start', color
    while l:
	eis,cc = l.pop(0)
	for e,i in eis:
	    if color[i] == S:
		if e != beta or vertex_color(cc[0])[i] != I:
		    #print ' reject'
		    memo[v] = None
		    return None
		color[i] = I
		#print ' ', color
	    elif color[i] == I:
		if e != gamma:
		    #print ' reject'
		    memo[v] = None
		    return None
		color[i] = S
		#print ' ', color
	    else:
		raise ValueError, 'Bad transition ' + str(eis) + ' in vertex_color'
    memo[v] = color
    #print 'is', color
    return color

def vertex_filter( v ):
    return vertex_color(v) is not None

def transition_filter( *tt ):
    return vertex_filter( tt[-1] )
