#*****************************************************************************
#  Copyright (C) 2017 Lee Worden <worden dot lee at gmail dot com>
#
#  Distributed under the terms of the GNU General Public License (GPL) v.2
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.all import *
from dynamicalsystems import *
from boxmodel import *
from boxmodelproduct import *

## inclusions:
##  given compartment 'c' of the i'th factor model,
##  and a compartment 'tup' of the product model,
##  is tup involved in the transitions of c in the product model, and
##  on which level(s).
## practically, will tup be allowed to take the place
##  of c in transitions.
## free_inclusions lets all free model compartments play all roles,
##  though it doesn't cross levels like the pair model does.
def free_inclusions( c, tup, i ):
    return [i]

def free_uop( s_tuple, iset, eis, transition_filter=lambda *x:True, edge_names={} ):
    def edge_name( r, exclude=() ):
        if r in edge_names: return edge_names[r]
        pp = [ p for p in r.variables() if p not in exclude ]
        return pp[0]
    s_eis = tuple( (edge_name(r,(v,w)),i) for (v,w,r),i in eis )
    #print 'uop', s_tuple, '+', (s_eis,())
    n_tuple = s_tuple + ( (s_eis,()), )
    if transition_filter( s_tuple, iset, eis, n_tuple ):
	#print ': ', n_tuple
	return [ n_tuple ]
    return []

def free_bop( s_tuple, iset, c_tuple, i_set, eis, transition_filter=lambda *x:True, edge_names={} ):
    def edge_name( r, exclude=() ):
        if r in edge_names: return edge_names[r]
        pp = [ p for p in r.variables() if p not in exclude ]
        return pp[0]
    s_eis = tuple( (edge_name(r,(v,w)),i) for (v,w,r),i in eis )
    #print 'bop', s_tuple, '+', (s_eis,(c_tuple,)) 
    n_tuple = s_tuple + ( (s_eis,(c_tuple,)), )
    if ( all( e == beta for e,i_ in s_eis ) and
         transition_filter( s_tuple, iset, c_tuple, i_set, s_eis, n_tuple )
       ):
	#print ': (', n_tuple, ')'
	return [ n_tuple ]
    return []

fw_memo_names={} 
def free_wrapper( *xx ):
    try: v = fw_memo_names[xx]
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
		        [ latex( subscriptedsymbol( eis[0][0], *(i for _,i in eis) ) ) ] +
		        [ '('+latex(free_wrapper(*c))+')'
                            if len(c) > 1 else
                            ' '+latex(free_wrapper(*c))
                            for c in cs
                        ]
                    )
                    for eis,cs in xx[1:]
	        ]
	    )
	    v = SR.symbol()
	    v = SR.symbol( str(v), latex_name=lname )
	    #print 'is', lname
	fw_memo_names[xx] = v
    return v

def free_param_relabel( *xx ):
    ## p, V, iota, W
    ## p, V, iota, C, iota_, W
    ## p, V, iota, iota_, W
    if len(xx) == 4:
        ## unary interaction with source compartment
        p, V, iota, W = xx
        return bm_param( p, *(iota + [free_wrapper(*V)]) )
    elif len(xx) == 5:
        ## binary interaction with interaction within source compartment
        p, V, iota, iota_, W = xx
        return bm_param( p, *(iota + [free_wrapper(*V)]) )
    elif len(xx) == 6:
        ## binary interaction with source interacting with catalyst
        p, V, iota, C, iota_, W = xx
        return bm_param( p, *(iota + [free_wrapper(*V), free_wrapper(*C)]) )

p_names = {}
def free_param_namer( *xx ):
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
def free_pos( graph, models, raw_tuples, init_pos={} ):
    ## how to position these compartments.
    ## initial boxes are positioned somewhat smartly
    ## every transition in level 0 adds 1 to x, and
    ## every transition in level 1 adds -1 to y.
    ## overlapping locations are shifted by (0.1,0.1) until unique.
    x_init_pos = {}
    angle = { 0:(0,-1), 1:(1,0), 2:(sqrt(3)/2,1/2), 3:(0,-3) }
    for i,m in enumerate(models):
        for v,mp in m._graph.get_pos().iteritems():
            x_init_pos[ (v,) ] = tuple( [ mx*ax for mx,ax in zip(mp,angle[i]) ] )
    x_init_pos.update( init_pos )
    init_pos = x_init_pos
    print 'init_pos', init_pos
    def free_pos_compartment( xx ):
	try: p = memo_pos[xx]
	except KeyError:
	    print 'pos', xx
	    if len(xx) == 1:
		#p = { S:(0,0), I_0:(1,0), I_1:(0,-1), I_01:(1,-1), I_10:(1,-1) }[ xx[0] ]
                p = init_pos[ (xx[0],) ]
		side = p[0] + p[1]
	    else:
		yx, eiscs = xx[:-1], xx[-1]
		free_pos_compartment( yx )
		p = vector( memo_pos_base[yx] )
		side = p[1] + p[0]
		print ' from', yx, 'at', p
		for _,i in eiscs[0]:
                    step = vector( ( (1,0), (0,-1) )[i] )
		    p = p + step
		    if p[0] + p[1] != 0:
                        side = p[0] + p[1]
		    print ' shift to', p
                p = tuple(p)
	    print ' base', p
	    memo_pos_base[xx] = p
	    while p in memo_pos_places:
		if side > 0:
	            p = ( p[0] + 0.1, p[1] + 0.15 )
		else:
	            p = ( p[0] - 0.1, p[1] - 0.15 )
	    print 'is', p
	    memo_pos[xx] = p
	    memo_pos_places[ p ] = 1
	return p
    return { free_wrapper(*xx) : free_pos_compartment(xx) for xx in raw_tuples }

def limited_uop( st, transition_filter=lambda *x:True ):
    def llu( s_tuple, iset, eis ):
        if s_tuple in st: return free_uop( s_tuple, iset, eis, transition_filter=transition_filter )
        else: return []
    return llu

def limited_bop( st, transition_filter=lambda *x:True ):
    def llb( s_tuple, iset, c_tuple, i_set, eis ):
        if s_tuple in st and c_tuple in st: return free_bop( s_tuple, iset, c_tuple, i_set, eis, transition_filter=transition_filter )
        else: return []
    return llb

## for use as a transition filter
def limit_to_iteration( n ):
    iteration_dict = {}
    def lim_it( *xx ):
        ## s_tuple, iset, c_tuple, i_set, s_eis, t_tuple
        s_tuple = xx[0]
        t_tuple = xx[-1]
        try:
            s_it = iteration_dict[s_tuple]
        except KeyError:
            s_it = 0
            iteration_dict[s_tuple] = s_it
        t_it = s_it + 1
        iteration_dict[t_tuple] = t_it
        return ( t_it <= n )
    return lim_it

def free_product( *models, **kwargs ):
    inclusions =        kwargs.pop( 'inclusions', free_inclusions )
    transition_filter = kwargs.pop( 'transition_filter', None )
    max_iteration =     kwargs.pop( 'max_iteration', None )
    if transition_filter is None and max_iteration is not None:
        transition_filter = limit_to_iteration( max_iteration )
    unary_operation =   kwargs.pop( 'unary_operation', lambda *x: free_uop(*x, transition_filter=transition_filter) )
    binary_operation =  kwargs.pop( 'binary_operation', lambda *x: free_bop(*x, transition_filter=transition_filter) )
    seed_set =          kwargs.pop( 'seed_set', None )
    compartment_renaming = kwargs.pop( 'compartment_renaming', lambda *x: x )
    compartment_wrapper = kwargs.pop( 'compartment_wrapper', free_wrapper )
    vertex_namer =      kwargs.pop( 'vertex_namer', lambda x:x )
    param_relabeling =  kwargs.pop( 'param_relabeling', free_param_relabel )
    param_namer =       kwargs.pop( 'param_namer', free_param_namer )
    vertex_positioner = kwargs.pop( 'vertex_positioner', free_pos )
    if seed_set is not None:
        def mktuple( s ):
            try:
                s[0]
            except TypeError:
                s = (s,)
            return s
        seed_set = [ mktuple(s) for s in seed_set ]
    return strong_product( *models,
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
	single_edge_generator = lambda eis, *a, **g: default_strong_edge_bundle_generator( eis, *a, **g ) if len( Set( (r for (v,w,r),i in eis) ) ) == 1 else []
    )

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
