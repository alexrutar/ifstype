from sympy import simplify, PurePoly, Rational as sy_R
from sympy.polys.domains import QQ
from sympy.polys.numberfields import minpoly
from sympy.polys.polytools import gcdex
from sympy.abc import x

from fractions import Fraction
import functools

#  from .interval import Interval, NetInterval

# -------------------------------------
# internal Rational class
# -------------------------------------
Rational = Fraction

# -------------------------------------
# constants
# -------------------------------------
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
    #  @constant
    #  def net_iv_base(self):
        #  return NetInterval(self.n_0,self.n_1,self.n_base)
    #  @constant
    #  def iv_0_1(self):
        #  return Interval(self.n_0,self.n_1)

Constants = _Const()



def cnst_poly(q): # workaround for sympy bug
    return PurePoly(x+q,domain=QQ)-PurePoly(x,domain=QQ)

def check_eq_anf(f):
    "Decorator to verify that algebraic numbers originate from equivalent factories, and convertes other to correct format."
    def valid_anf(self,other):
        if isinstance(other,Rational):
            return f(self, cnst_poly(sy_R(other.numerator, other.denominator)))
        elif isinstance(other,int): # for powers
            return f(self, other)
        elif isinstance(other,AlgebraicNumber) and self.anf == other.anf:
            return f(self, other._poly)
        else:
            raise ValueError("Inequivalent algebraic numbers!")
    return valid_anf

# TODO: write this later
class RationalModPoly:
    """Rational polynomial class which supports multiplication and inversion on tuples, modulo a fixed irreducible polynomial.
    Internally, polynomials are just tuples, (a0,a1,...,an) ~ a0+a1*x+...+an*x^n"""
    def __init__(self, minpoly):
        self.poly = poly
        self.deg = len(poly)

# algebraic numbers
class AlgebraicNumber:
    def __init__(self, anf, poly):
        self.anf = anf # pointer to the factory
        self._poly = poly # sympy PurePoly

    def __repr__(self):
        return f"x={self.anf.expr}, expr={self._poly}"
    def __str__(self):
        return f"{str(self._poly.as_expr())}"

    def __hash__(self):
        return hash(tuple(self._poly.coeffs()))

    # -----------------------------------------
    # comparisons
    # -----------------------------------------
    @check_eq_anf
    def __lt__(self, other):
        return self.anf.lt(self,other)

    @check_eq_anf
    def __le__(self, other):
        return self.anf.le(self,other)

    @check_eq_anf
    def __gt__(self, other):
        return self.anf.gt(self,other)

    @check_eq_anf
    def __ge__(self, other):
        return self.anf.ge(self,other)

    @check_eq_anf
    def __eq__(self,other):
        return self.anf.eq(self,other)

    @check_eq_anf
    def __neq__(self,other):
        return self.anf.neq(self,other)
    # -----------------------------------------
    # numeric
    # -----------------------------------------

    @check_eq_anf
    def __truediv__(self, other):
        return self.anf.mul(self,self.anf.rinv(other)._poly)

    @check_eq_anf
    def __rtruediv__(self, other):
        return self.anf.mul(self.anf.inv(self),other)

    @check_eq_anf
    def __add__(self,other):
        return self.anf.add(self,other)

    @check_eq_anf
    def __radd__(self,other):
        return self.anf.add(self,other)

    @check_eq_anf
    def __sub__(self,other):
        return self.anf.sub(self,other,flip=False)

    @check_eq_anf
    def __rsub__(self,other):
        return self.anf.sub(self,other,flip=True)

    @check_eq_anf
    def __mul__(self,other):
        if type(other) is float:
            print(other)
        return self.anf.mul(self,other)

    @check_eq_anf
    def __rmul__(self,other):
        return self.anf.mul(self,other)
        
    def __pow__(self, it):
        return self.anf.pow(self,it)

    # ---------------------------------
    # unary operations
    # ---------------------------------
    def __abs__(self):
        return self.anf.abs(self)
    def __neg__(self):
        return self.anf.neg(self)
    def __pos__(self):
        return self.anf.post(self)

    def __float__(self):
        return self.anf.float(self)

class AlgebraicNumberFactory:
    """A class governing the creation of algebraic numbers.
    Two factories are considered equivalent if they come from the same expression.
    Represents algebraic numbers with respect to their basis over Q."""
    def __init__(self, expr):
        self.expr = simplify(expr)
        self.min_poly = minpoly(expr,polys=True,domain=QQ)
        self.deg = self.min_poly.degree()

    #  @functools.lru_cache
    def simplify_pow(self,alpha:int):
        pass

    def __eq__(self, other):
        return self.expr == other.expr

    def _from_poly(self, poly):
        return AlgebraicNumber(self, poly.rem(self.min_poly))

    def one(self):
        return self._from_poly(cnst_poly(1))

    def alpha(self):
        return self._from_poly(PurePoly(x,domain=QQ))

    def inv(self, an):
        return self._from_poly(an._poly.invert(self.min_poly))
    def rinv(self, q_or_poly):
        if isinstance(q_or_poly,PurePoly):
            out = self._from_poly(q_or_poly.invert(self.min_poly))
            return out
        else:
            return q_or_poly**(-1)

    # -----------------------------------------
    # comparisons
    # -----------------------------------------
    def lt(self, an, q_or_poly):
        if isinstance(q_or_poly,PurePoly):
            q_or_poly = q_or_poly.eval(self.expr)
        return an._poly.eval(self.expr) < q_or_poly

    def le(self, an, q_or_poly):
        if isinstance(q_or_poly,PurePoly):
            q_or_poly = q_or_poly.eval(self.expr)
        return an._poly.eval(self.expr) <= q_or_poly

    def gt(self, an, q_or_poly):
        if isinstance(q_or_poly,PurePoly):
            q_or_poly = q_or_poly.eval(self.expr)
        return an._poly.eval(self.expr) > q_or_poly

    def ge(self, an, q_or_poly):
        if isinstance(q_or_poly,PurePoly):
            q_or_poly = q_or_poly.eval(self.expr)
        return an._poly.eval(self.expr) >= q_or_poly

    def eq(self, an, q_or_poly):
        if isinstance(q_or_poly,Rational):
            q_or_poly = cnst_poly(q_or_poly)
        return an._poly == q_or_poly

    def neq(self, an, q_or_poly):
        if isinstance(q_or_poly,Rational):
            q_or_poly = cnst_poly(q_or_poly)
        return an._poly != q_or_poly
    # -------------------
    # __add__, __radd__
    # -------------------
    def add(self, an1, q_or_poly):
        return self._from_poly(an1._poly.add(q_or_poly))

    # -------------------
    # __mul__, __rmul__
    # -------------------
    def mul(self, an1, q_or_poly):
        #  print(f"an1={an1},q_or_poly={q_or_poly}")
        #  print("normal prod ",an1._poly*q_or_poly)
        #  print("internal mul ",an1._poly.mul(q_or_poly))
        return self._from_poly(an1._poly.mul(q_or_poly))

    # -------------------
    # __sub__, __rsub__
    # -------------------
    def sub(self, an1, q_or_poly,flip=False):
        if flip:
            return self._from_poly(-an1._poly.sub(q_or_poly))
        else:
            return self._from_poly(an1._poly.sub(q_or_poly))

    def pow(self, an, it):
        if it < -1:
            raise NotImplementedError("Power must be integer >= -1")
        elif it == -1:
            return self.inv(an)
        elif it == 0:
            return self.one()
        elif it == 1:
            return self._from_poly(an._poly)
        elif it == 2:
            return self.mul(an,an._poly)
        # use divide and conquer algorithm for larger powers
        else:
            return self.mul(self.pow(an,-(-it // 2)),self.pow(an,(it // 2))._poly)

    # -------------------
    # unary operations
    # -------------------
    def abs(self, an):
        val = an._poly.eval(self.expr)
        if val < 0:
            return self._from_poly(-an._poly)
        else:
            return self._from_poly(an._poly)

    def neg(self, an):
        return self._from_poly(-an._poly)
    def pos(self, an):
        return self._from_poly(an._poly)

    # -------------------
    # __float__
    # -------------------
    def float(self, an):
        return float(an._poly.eval(self.expr))

    
