#*****************************************************************************
#  Copyright (C) 2017 Lee Worden <worden dot lee at gmail dot com>
#
#  Distributed under the terms of the GNU General Public License (GPL) v.2
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#from sage.all import *
from boxmodelproduct import *
from dynamicalsystems import *

## code evaluating categories of compartments in R

def R_inclusions_fn( self, name='observations', inclusions=None, extras=Bindings() ):
    """R_inclusions_fn: emit definition of an R function that constructs
    aggregate quantities from the compartments of a product model.

    inclusions: which quantities to define, if not the ones generated in
    the process of the product operation by tracking the division of
    factor compartments into product compartments.
    extras: quantities to include in addition to the above."""
    code = '#!/usr/bin/R\n'
    code += name + ' <- function( state ) {\n'
    code += '  with(state, {\n'
    code += '    obs <- list(c(\n'
    if inclusions is None:
        code += ',\n'.join(
            '      ' + str(v) + ' = ' + ' + '.join( str(vt) for vt in ll )
            for v, ll in self._inclusion_variables.iteritems()
            if ll != [v]
        ) + '\n'
    else:
        code += ',\n'.join(
            '      ' + str(k) + ' = ' + str(v)
            for k, v in inclusions._dict.iteritems()
        ) + '\n'
    if len(extras._dict) > 0:
        code += ',\n'.join(
            '      ' + str(k) + ' = ' + str(v)
            for k, v in extras._dict.iteritems()
        ) + '\n'
    code += '    ))\n'
    code += '  })\n'
    code += '  return(obs)\n'
    code += '}\n'
    return code

BoxModelProduct.R_inclusions_fn = R_inclusions_fn

def R_marginal_names( self, name='marginals' ):
    """R_marginal_names_fn: provide for R the names of compartments
    indexed by compartments of the factor models"""
    code = '#!/usr/bin/R\n'
    code += name + ' <- c(\n'
    if len(self._variable_marginals):
        code += ',\n'.join(
            '  ' + str(v) + ' = c("' + '", "'.join( str(vt) for vt in ll ) + '")'
            for v, ll in self._variable_marginals.iteritems()
            if ll != [v]
        ) + '\n'
    if len(self._parameter_marginals):
        code += ',\n'.join(
            '  ' + str(v) + ' = c("' + '", "'.join( str(vt) for vt in ll ) + '")'
            for v, ll in self._parameter_marginals.iteritems()
            if ll != [v]
        ) + '\n'
    code += ')\n'
    return code

BoxModelProduct.R_marginal_names = R_marginal_names
