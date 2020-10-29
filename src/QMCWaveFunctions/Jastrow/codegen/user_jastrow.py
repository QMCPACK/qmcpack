from __future__ import print_function

from sympy import Symbol, diff, exp, Piecewise
from sympy.printing.cxxcode import CXX11CodePrinter
import sys

# Generate code for a radial function for Jastrow factors from a symbolic expression.
# This script creates ../UserFunctor.h from UserFunctor.h.in.
# It also prints the expression, derivatives, and a sample XML input block.


# The radial function is a function of r
r = Symbol('r')


# Possible symbols for variational parameters or fixed values
A = Symbol('A')
B = Symbol('B')
C = Symbol('C')
D = Symbol('D')
rc = Symbol('R_c',positive=True)


# optimizable parameters
variational_parameters = []

# parameters from the input file, but are not to be optimized
input_parameters = []

# For now, assume the cusp value is one of the parameters (or None)
cusp_param = None

# Define the functional form below.

# Simple Pade
if True:
  f = A*r/ (1+B*r) - A/B
  variational_parameters = [A,B]
  cusp_param = A

# second order Pade
if False:
  f = (A*r + C*r*r)/(1 + B*B*r + D*D*r*r)
  variational_parameters = [A,B,C,D]

# Gaussian
if False:
  f1 = A*exp(-r**2/C**2)
  # For smooth truncation
  f2 = f1.subs(r, 2*rc - r) - 2*f1.subs(r, rc)
  f = Piecewise((f2, r < rc), (0, True))

  variational_parameters = [A,C]
  input_parameters = [rc]

# Wagner-Mitas  https://arxiv.org/abs/cond-mat/0610088
if False:
  x = Symbol('x')
  z =  x**2*(6 - 8*x + 3*x*x)
  f = (1 - z.subs(x, r/rc))/(1 + B*z.subs(x, r/rc))
  variational_parameters = [B]
  input_parameters = [rc]

# McMillan Phys Rev 138 442
if False:
  f = (B/r)**A
  variational_parameters = [A,B]


#------------------------
# End of user input area
#------------------------


def check_free_symbols(f, variational_parameters, input_parameters):
  # Check all free symbols are accounted for in the parameter lists
  free_sym = f.free_symbols
  all_sym = set(variational_parameters + input_parameters + [r])
  missing_in_params = free_sym - all_sym
  if missing_in_params:
    print('Symbols are missing from variational or input parameter lists: ',missing_in_params)
    sys.exit(1)

  extra_in_params = all_sym - free_sym
  if extra_in_params:
    print('Extra symbols in variational or input parameter lists: ',extra_in_params)
    sys.exit(1)






# Convert squares to multiplications

class NoPowCodePrinter(CXX11CodePrinter):
  def __init__(self, settings=None):
    super(NoPowCodePrinter, self).__init__(settings=settings)

  def _print_Pow(self, expr):
    if expr.exp == 2:
      e = self._print(expr.base)
      return '((%s)*(%s))'%(e,e)
    return super(NoPowCodePrinter, self)._print_Pow(expr)

CP = NoPowCodePrinter()


def get_var_names(var_name):
  """Defines replacement dictionary for the bare variable name and
     the names derived from it - the optimization flag and the identifier name.
  """
  repl = dict()
  repl['opt_var_name'] = "Opt_%s"%var_name
  repl['id_var_name'] = "ID_%s"%var_name
  repl['var_name'] = var_name
  return repl

# String representation of the function, for the comments
def gen_func_str(f):
  return ' *  f = %s\n'%f


# Parameter declarations in the class
def gen_param_defs(param_list, input_param=None):
  out_str = """
  /// Is optimizable
  bool %(opt_var_name)s;
  /// Value
  real_type %(var_name)s;
  /// Id in XML input
  std::string %(id_var_name)s;
  """

  s = ""
  for var_name in param_list:
    #opt_var = "Opt_%s"%var_name
    #id_var = "ID_%s"%var_name
    repl = get_var_names(var_name)
    s += out_str%repl
    s += "\n"

  if input_param:
    for var_name in input_param:
      repl = get_var_names(var_name)
      s += out_str%repl
      s += "\n"

  print("Creating parameter variable definitions")
  print(s)
  print("")

  return s


# Function for setting the cusp value
def gen_set_cusp(cusp_param):
  out_str = """
  void setCusp(real_type cusp)
  {
    %(var_name)s     = cusp;
    %(opt_var_name)s = false;
    reset();
  }
  """

  s = ""
  if cusp_param:
    repl = get_var_names(cusp_param)
    s = out_str%repl
    print("Set cusp function")
    print(s)
  else:
    print("No cusp")

  print()

  return s


# Evaluate function only
def gen_evaluate(f):
  out_str = """
  inline real_type evaluate(real_type r) const {
   return %s;
  }
  """

  val_str = CP.doprint(f)
  s = out_str%val_str
  print('Creating evaluate function')
  print(s)
  print('')
  return s


# Evaluate function plus first and second derivatives
def gen_evaluate_2nd_deriv(f, df, ddf):

  out_str = """
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2) const
  {
    dudr        = %(df)s;
    d2udr2      = %(ddf)s;
    return %(val)s;
  }

"""
  val_str = CP.doprint(f)
  df_str = CP.doprint(df)
  ddf_str = CP.doprint(ddf)
  s = out_str%{"val":val_str,  "df":df_str, "ddf":ddf_str}
  print('Creating evaluate 1st and 2nd derivatives function')
  print(s)
  print('')

  return s


# Evaluate function plus first, second, and third derivatives
def gen_evaluate_3rd_deriv(f, df, ddf, d3f):

  out_str = """
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type& d3udr3) const
  {
    dudr        = %(df)s;
    d2udr2      = %(ddf)s;
    d3udr3      = %(d3f)s;
    return %(val)s;
  }

"""
  val_str = CP.doprint(f)
  df_str = CP.doprint(df)
  ddf_str = CP.doprint(ddf)
  d3f_str = CP.doprint(d3f)
  s = out_str%{"val":val_str,  "df":df_str, "ddf":ddf_str, "d3f":d3f_str}
  print('Creating evaluate 1st, 2nd, and 3rd derivatives function')
  print(s)
  print('')

  return s


# Evaluate derivatives of just the function with respect to the parameters
def gen_evaluate_parameter_derivatives(variational_parameters, param_derivs):

  func_out_str = """
  /// compute derivatives with respect to variational parameters
  inline bool evaluateDerivatives(real_type r, std::vector<real_type>& derivs)
  {
    int i = 0;
    %s

    return true;
  }
  """

  each_block_str = """
    if (%(opt_var_name)s)
    {
      derivs[i] = %(param_deriv_value)s;
      ++i;
    }
  """

  block_s = ""
  for vp in variational_parameters:
    repl = get_var_names(vp)
    repl['param_deriv_value'] = CP.doprint(param_derivs[vp][0])
    block_s += each_block_str%repl

  s = func_out_str%block_s


  print('Creating evaluate parameter derivative function')
  print(s)
  print('')

  return s


# Evaluate derivatives of the function (and spatial derivatives) with respect to the parameters
def gen_evaluate_all_parameter_derivatives(variational_parameters, param_derivs):
  func_out_str = """
  inline bool evaluateDerivatives(real_type r, std::vector<TinyVector<real_type, 3>>& derivs)
  {
    int i = 0;
    %s
    return true;
  }
  """

  each_block_str = """
    if (%(opt_var_name)s)
    {
      derivs[i][0] = %(param_deriv_value)s;
      derivs[i][1] = %(param_deriv_first)s;
      derivs[i][2] = %(param_deriv_second)s;
      ++i;
    }
  """

  block_s = ""
  for vp in variational_parameters:
    repl = get_var_names(vp)
    repl['param_deriv_value'] = CP.doprint(param_derivs[vp][0])
    repl['param_deriv_first'] = CP.doprint(param_derivs[vp][1])
    repl['param_deriv_second'] = CP.doprint(param_derivs[vp][2])
    block_s += each_block_str%repl

  s = func_out_str%block_s


  print('Creating evaluate all parameter derivatives function')
  print(s)
  print('')

  return s

#  Handle the XML input
def gen_xml_input(param_list, input_param_list=None):
  out_str="""
  bool put(xmlNodePtr cur)
  {
    cur = cur->xmlChildrenNode;
    while (cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      if (cname == "var") //only accept var
      {
        std::string id_in;
        std::string p_name;
        OhmmsAttributeSet rAttrib;
        rAttrib.add(id_in, "id");
        rAttrib.add(p_name, "name");
        rAttrib.put(cur);
        %(read_param)s
      }
      cur = cur->next;
    }
    reset();
    myVars.clear();
    %(set_opt)s
    return true;
  }
"""

  read_param_str = """
        if (p_name == "%(var_name)s")
        {
          %(id_var_name)s = id_in;
          putContent(%(var_name)s, cur);
          %(opt_var_name)s = %(do_opt)s;
        }
  """

  set_opt_str = """
    if (%(opt_var_name)s)
      myVars.insert(%(id_var_name)s, %(var_name)s, %(opt_var_name)s, optimize::OTHER_P);
  """

  read_param = ""
  set_opt = ""
  for var_name in param_list:
    repl = get_var_names(var_name)
    repl['do_opt'] = 'true'
    read_param += read_param_str%repl
    set_opt += set_opt_str%repl

  if input_param_list:
    for var_name in input_param_list:
      repl = get_var_names(var_name)
      repl['do_opt'] = 'false'
      read_param += read_param_str%repl


  s = out_str%{"read_param":read_param, "set_opt":set_opt}

  print('Creating XML input function')
  print(s)
  print('')

  return s


# Code for resetting the parameters
def gen_reset_parameters(param_list):
  out_str = """
  void resetParameters(const opt_variables_type& active)
  {
    if (myVars.size())
    {
      int ia = myVars.where(0);
      if (ia > -1)
      {
        int i = 0;
        %(set_param)s
      }
      reset();
    }
  }
  """

  set_param = """
        if (%(opt_var_name)s)
          %(var_name)s = std::real(myVars[i++] = active[ia++]);
  """

  set_param_str = ""
  for var_name in param_list:
    repl = get_var_names(var_name)
    set_param_str += set_param%repl

  s = out_str%{"set_param":set_param_str}

  print('Creating reset parameters function')
  print(s)
  print('')

  return s


# Example XML input
def gen_example_xml_input(param_list, input_param_list=None, cusp_param=None):
  out_str = """
   <jastrow name="Jee" type="Two-Body" function="user">
      <correlation speciesA="u" speciesB="d">
%(var_entries)s
      </correlation>
    </jastrow>
   """
  one_entry_str = """        <var id="j_%(var_name)s" name="%(var_name)s">1.0</var>\n"""

  entries = ""

  for var_name in param_list:
    if var_name != cusp_param:
      entries += one_entry_str%{'var_name':var_name}


  if input_param_list:
    for var_name in input_param_list:
      entries += one_entry_str%{'var_name':var_name}

  s = out_str%{'var_entries':entries}

  print("Example XML Input section")
  print(s)
  print()



dire_codegen_text = """
/*
 DO NOT MAKE PERMANENT EDITS IN THIS FILE
 This file is generated from src/QMCWavefunctions/Jastrow/codegen/user_jastrow.py and UserFunctor.h.in

 To make changes, edit UserFunctor.h.in and rerun user_jastrow.py.
*/
"""

def run_template(fname_in, fname_out, bodies):
  from string import Template
  out = ''
  with open(fname_in, 'r') as f:
    s = Template(f.read())
    try:
      out = s.substitute(bodies)
    except KeyError as e:
      print('Error, template item not found: ',e)
    with open(fname_out, 'w') as fout:
      fout.write(out)

def run_all(f, variational_parameters, input_parameters, cusp_param):
  print("Input for radial Jastrow:")
  print(f)

  check_free_symbols(f, variational_parameters, input_parameters)

  # Compute radial derivatives

  df = diff(f, r)
  ddf = diff(f, r, 2)
  d3f = diff(f, r, 3)

  print()
  print("Cusp (df/dr at r = 0) = ",df.subs(r,0))
  print()

  param_derivs = dict()
  for vp in variational_parameters:
    pd_f = diff(f, vp)
    pd_df = diff(df, vp)
    pd_ddf = diff(ddf, vp)
    param_derivs[vp] = (pd_f, pd_df, pd_ddf)

  bodies = {
    'dire_codegen_warning' : dire_codegen_text,
    'func_str' : gen_func_str(f),
    'param_defs' : gen_param_defs(variational_parameters, input_parameters),
    'set_cusp' : gen_set_cusp(cusp_param),
    'evaluate_func' : gen_evaluate(f),
    'evaluate_func_2nd_derivative' : gen_evaluate_2nd_deriv(f, df, ddf),
    'evaluate_func_3rd_derivative' : gen_evaluate_3rd_deriv(f, df, ddf, d3f),
    'evaluate_parameter_derivative' : gen_evaluate_parameter_derivatives(variational_parameters, param_derivs),
    'evaluate_all_parameter_derivatives' : gen_evaluate_all_parameter_derivatives(variational_parameters, param_derivs),
    'xml_input' : gen_xml_input(variational_parameters, input_parameters),
    'reset_parameters' : gen_reset_parameters(variational_parameters)
  }

  run_template('UserFunctor.h.in', '../UserFunctor.h', bodies)

  gen_example_xml_input(variational_parameters, input_parameters, cusp_param)


if __name__ == '__main__':
  run_all(f, variational_parameters, input_parameters, cusp_param)
