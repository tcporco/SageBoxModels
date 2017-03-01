#from sage.all import *
from boxmodelproduct import *
from dynamicalsystems import *

## code evaluating categories of compartments in R

def R_inclusions_fn( self, name='observations', extras=Bindings() ):
    code = '#!/usr/bin/R\n'
    code += name + ' <- function( state ) {\n'
    code += '  with(state, {\n'
    code += '    obs <- list(c(\n'
    code += ',\n'.join(
            '      ' + str(v) + ' = ' + ' + '.join( str(vt) for vt in ll )
            for v, ll in self._inclusion_variables.iteritems()
            if ll != [v]
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
