from __future__ import print_function

# Extensions to Sympy code generation to support C++ with templates and more

# Used in gen_cubic_spline_solver.py

from sympy import Set, Basic, Tuple, IndexedBase
from sympy.codegen.ast import Assignment, Pointer, Node, Type
from sympy.codegen.ast import String, Declaration, Variable
from sympy.printing.cxxcode import CXX11CodePrinter

# A range class that accepts symbolic limits. Purpose is For loops
class ARange(Set):
    is_iterable = True

    def __new__(cls, *args):
        slc = slice(*args)
        start = slc.start
        stop = slc.stop
        step = slc.step

        return Basic.__new__(cls, start, stop, step)
    start = property(lambda self: self.args[0])
    stop = property(lambda self: self.args[1])
    step = property(lambda self: self.args[2])

    # Just here to make it pass the 'iterable' test
    def __iter__(self):
        i = 0
        yield i

# The end value should be included in the iteration
class ARangeClosedEnd(ARange):
  pass


# Convert Eq to Assignment
def convert_eq_to_assignment(expr):
    return Assignment(expr.lhs, expr.rhs)


# Node for a C++ reference
class Reference(Pointer):
    """ Represents a C++ reference"""
    pass


# Specify direct initialization (with parentheses)
#  e.g.   int j(0);
class VariableWithInit(Variable):

    __slots__ = ['type_init'] + Variable.__slots__


# Templated function definition
class TemplateFunctionDefinition(Node):
    __slots__ = ['return_type','name','parameters','template_types','body','attrs']
    _construct_return_type = Type
    _construct_name = String

    @staticmethod
    def _construct_parameters(args):
        def _var(arg):
            if isinstance(arg, Declaration):
                return arg.variable
            elif isinstance(arg, Variable):
                return arg
            else:
                return Variable.deduced(arg)
        return Tuple(*map(_var, args))

    @staticmethod
    def _construct_template_types(args):
        return Tuple(*args)


# Code printer for extended features
class ACodePrinter(CXX11CodePrinter):
    def __init__(self, settings=None):
        super(ACodePrinter, self).__init__(settings=settings)

    def _print_Assignment(self, expr):
        lhs = expr.lhs
        rhs = expr.rhs
        if lhs.has(IndexedBase) or rhs.has(IndexedBase):
            return self._get_statement("%s = %s"%(self._print(lhs),self._print(rhs)))
        else:
            return super(ACodePrinter, self)._print_Assignment(expr)

    def _print_Pow(self, expr):
      if expr.exp == 2:
        e = self._print(expr.base)
        return '%s*%s'%(e,e)
      return super(ACodePrinter, self)._print_Pow(expr)

    def _print_Symbol(self, expr):
        name = super(ACodePrinter, self)._print_Symbol(expr)
        # Replace prime marker in symbol name with something acceptable in C++
        #  Maybe should generalize to a lookup from symbol name to code name?
        if "'" in name:
            name = name.replace("'", "p")
        return name

    def _print_Declaration(self, decl):
        #print("decl = ",decl,type(decl))
        var = decl.variable
        val = var.value
        if isinstance(var, Reference):
            result = '{t}& {s}'.format(
                t = self._print(var.type),
                s = self._print(var.symbol)
            )
            return result
        elif isinstance(var, VariableWithInit):
            result = '{t} {s}({init})'.format(
                      t=self._print(var.type),
                      s=self._print(var.symbol),
                      init=self._print(var.type_init))
            return result
        else:
            return super(ACodePrinter, self)._print_Declaration(decl)

    def _print_TemplateFunctionDefinition(self, expr):
        decl = "template<{template_args}>\n{ret_type} {name}({params}){body}".format(
                    template_args=', '.join(map(lambda arg: 'typename '+self._print(arg), expr.template_types)),
                    ret_type=self._print(expr.return_type),
                    name=expr.name,
                    params=', '.join(map(lambda arg: self._print(Declaration(arg)), expr.parameters)),
                    body=self._print_Scope(expr)

                )
        return decl

    def _print_For(self, expr):
        target = self._print(expr.target)
        it = expr.iterable
        body = self._print(expr.body)
        #print("it = ",it,type(it),isinstance(it,ARange))

        if isinstance(it, ARange):
            end_compare = ""
            if isinstance(it, ARangeClosedEnd):
              end_compare="="
            if it.step > 0:
                return ("for (auto {target} = {start}; {target} <{end_compare} {stop}; {target} += {step}) {{\n{body}\n}}").format(
                    target=target,start=it.start, stop=it.stop, step=it.step, body=body,end_compare=end_compare)
            else:
                return ("for (auto {target} = {start}; {target} >{end_compare} {stop}; {target} += {step}) {{\n{body}\n}}").format(
                    target=target,start=it.start, stop=it.stop, step=it.step, body=body, end_compare=end_compare)
        else:
            return super(ACodePrinter, self)._print_For(expr)






