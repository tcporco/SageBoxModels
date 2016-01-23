import boxmodel

# in the future: create a StochasticProcess class, put most of this
# there, and have BoxModel produce one?

def R_ode_fn( self, name='odefn' ):
    vars = sorted(self._vars, key=str)
    params = sorted(self._parameters, key=str)
    code = '#!/usr/bin/R\n'
    code += name + ' <- function( t, x, pa ) {\n'
    code += '  ## state\n'
    for v in vars:
	code += '  ' + str(v) + " <- x[['" + str(v) + "']]\n"
    code += '  ## parameters\n'
    for p in params:
	code += '  ' + str(p) + " <- pa[['" + str(p) + "']]\n"
    code += '  ## the derivative\n'
    flow = self.ode_flow()
    code += '  dxdt <- list(c(\n'
    code += ',\n'.join( '    ' + str(v) + ' = ' + str(flow[v]) for v in vars ) + '\n'
    code += '  ))\n'
    code += '  print(dxdt)\n'
    code += '  return(dxdt)\n'
    code += '}\n'
    return code

boxmodel.BoxModel.R_ode_fn = R_ode_fn
