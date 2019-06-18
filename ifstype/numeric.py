from sympy import Rational, simplify
from sympy.polys.domains import QQ
from sympy.polys.numberfields import minpoly

from .interval import Interval, NetInterval

def as_relation(expr):
    poly_lst = minpoly(expr,polys=True,domain=QQ).all_coeffs()
    leading = poly_lst[0]
    return [p/leading for p in poly_lst]

# constants
def constant(f):
    def fset(self, value):
        raise AttributeError("Cannot change constant values")
    def fget(self):
        return f(self)
    return property(fget, fset)

class _Const:
    @constant
    def n_base(self):
        return Rational(2)
    @constant
    def n_0(self):
        return Rational(0)
    @constant
    def n_1(self):
        return Rational(1)
    @constant
    def net_iv_base(self):
        return NetInterval(self.n_0,self.n_1,self.n_base)
    @constant
    def iv_0_1(self):
        return Interval(self.n_0,self.n_1)

Constants = _Const()


class AlgebraicNumberFactory:
    """A class governing the creation of algebraic numbers.
    Two factories are considered equivalent if they come from the same expression.
    Represents algebraic numbers with respect to their basis over Q."""
    def __init__(self, expr):
        self.expr = simplify(expr)
        self.poly = minpoly(expr,polys=True,domain=QQ)
        self.deg = len(self.poly)-1

    def __eq__(self, other):
        return self.expr == other.expr

    def one(self):
        return self.from_basis(tuple(1)+tuple(Rational(0) for _ in range(self.deg-1)))

    def _from_basis(self, base):
        # internal method does not check length: base must have length self.deg
        # can take any generator expression
        return AlgebraicNumber(self, tuple(base))

    def add(self, an1, an2):
        return self._from_basis(a+b for a,b in zip(an1.basis, an2.basis))

    def sub(self, an1, an2):
        return self._from_basis(a-b for a,b in zip(an1.basis, an2.basis))

    #  def inv(self, an1)

def check_eq_anf(f):
    "Decorator to verify that algebraic numbers originate from equivalent factories"
    def valid_anf(self,other):
        if self.anf != other.anf:
            raise ValueError("Inequivalent algebraic numbers!")
        return f(self,other)
    return valid_anf
    
# algebraic numbers
class AlgebraicNumber:
    def __init__(self, anf, basis):
        self.anf = anf # pointer to the factory
        self.basis = basis # tuple

    def __str__(self):
        return f"val={self.anf.expr}, base={self.basis}"

    @check_eq_anf
    def __add__(self,other):
        return self.anf.add(self,other)

    @check_eq_anf
    def __sub__(self,other):
        return self.anf.sub(self,other)
