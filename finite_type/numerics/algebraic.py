import math
import itertools
import functools

from .polynomial import Poly
from .rational import Rational, Constants as C

# parse the args to convert a Rational argument to AlgebraicNumber, to be used with AlgebraicNumberFactory
def check_other(f):
    "Decorator to verify that algebraic numbers originate from equivalent factories, and convertes Rationals to correct format."
    def valid_anf(self,other):
        if isinstance(other,(int,Rational)):
            return f(self, self.anf._from_poly(Poly.cnst(other)))
        elif isinstance(other,AlgebraicNumber) and self.anf == other.anf:
            return f(self, other)
        else:
            raise ValueError(f"Inequivalent algebraic numbers! Other is a {type(other)}.")
    return valid_anf

class AlgebraicNumber:
    def __init__(self, anf, poly):
        self.anf = anf # pointer to the factory
        self._poly = poly

    def __repr__(self):
        return f"x={self.anf.expr_float}, expr={self._poly}"
    def __str__(self):
        return f"{self._poly.with_symbol('a')}"
    def as_latex(self,symb='a'):
        return f"{self._poly.with_symbol(symb)}"

    def __hash__(self):
        return hash(self._poly)

    # -----------------------------------------
    # comparisons
    # -----------------------------------------
    @check_other
    def __lt__(self, other):
        return self.anf.lt(self,other)

    @check_other
    def __le__(self, other):
        return self.anf.le(self,other)

    @check_other
    def __gt__(self, other):
        return self.anf.gt(self,other)

    @check_other
    def __ge__(self, other):
        return self.anf.ge(self,other)

    @check_other
    def __eq__(self,other):
        return self.anf.eq(self,other)

    @check_other
    def __neq__(self,other):
        return self.anf.neq(self,other)
    # -----------------------------------------
    # numeric
    # -----------------------------------------

    @check_other
    def __truediv__(self, other):
        return self.anf.mul(self,self.anf.inv(other))

    @check_other
    def __rtruediv__(self, other):
        return self.anf.mul(self.anf.inv(self),other)

    @check_other
    def __add__(self,other):
        return self.anf.add(self,other)

    @check_other
    def __radd__(self,other):
        return self.anf.add(self,other)

    @check_other
    def __sub__(self,other):
        return self.anf.sub(self,other)

    @check_other
    def __rsub__(self,other):
        return self.anf.sub(other,self)

    @check_other
    def __mul__(self,other):
        return self.anf.mul(self,other)

    @check_other
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
        return self.anf.pos(self)

    def __float__(self):
        return self.anf.float(self)


class AlgebraicNumberFactory:
    """A class governing the creation of algebraic numbers.
    Two factories are considered equivalent if they come from the same expression.
    Represents algebraic numbers with respect to their basis over Q.
    Internally, for comparisons, we use a float form for the expression. Be careful, if you need a lot of precision, this could be dangerous!"""
    # TODO: have arbitrary precision float / detect when floats may not be precise enough.
    # for most applications, this should be fine, at least in finite type IFS since all the comparisons will be with intervals about the same length
    def __init__(self, minpoly, expr_float):
        self.minpoly = minpoly
        self.expr_float = expr_float
        self.deg = self.minpoly.deg

    @classmethod
    def from_sympy(cls, sy_expr):
        from sympy.polys.domains import QQ
        from sympy.polys.numberfields import minpoly

        minpoly = Poly.from_PurePoly(minpoly(sy_expr,polys=True,domain=QQ))
        expr_float = float(sy_expr)
        return cls(minpoly,expr_float)

    def __eq__(self, other):
        return math.isclose(self.expr_float, other.expr_float)

    def _from_poly(self, poly):
        return AlgebraicNumber(self, poly.rem(self.minpoly))

    def zero(self):
        return self._from_poly(Poly.zero())

    def one(self):
        return self._from_poly(Poly.one())

    def alpha(self):
        return self._from_poly(Poly.monomial(1))

    def inv(self, an):
        return self._from_poly(an._poly.invert(self.minpoly))

    # -----------------------------------------
    # comparisons
    # exact comparison is precise, but checking for indirect comparison depends on the float value
    # -----------------------------------------
    def lt(self, an1, an2):
        return an1._poly.eval(self.expr_float) < an2._poly.eval(self.expr_float)

    def eq(self, an1, an2):
        return an1._poly == an2._poly

    def le(self, an1, an2):
        return self.eq(an1, an2) or self.lt(an1, an2)

    def gt(self, an1, an2):
        return self.lt(an2,an1)

    def ge(self, an1, an2):
        return self.eq(an1,an2) or self.lt(an2,an1)


    def neq(self, an1, an2):
        return an1._poly != an2._poly

    # -------------------
    # __add__, __radd__
    # -------------------
    def add(self, an1, an2):
        return self._from_poly(an1._poly + an2._poly)

    # -------------------
    # __mul__, __rmul__
    # -------------------
    def mul(self, an1, an2):
        return self._from_poly(an1._poly * an2._poly)

    # -------------------
    # __sub__, __rsub__
    # -------------------
    def sub(self, an1, an2):
        return self._from_poly(an1._poly - an2._poly)

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
            return self.mul(an,an)
        # use divide and conquer algorithm for larger powers
        else:
            return self.mul(self.pow(an,-(-it // 2)),self.pow(an,(it // 2)))

    # -------------------
    # unary operations
    # -------------------
    def abs(self, an):
        val = an._poly.eval(self.expr_float)
        if val < 0:
            return self.neg(an)
        else:
            return self.pos(an)

    def neg(self, an):
        return self._from_poly(-an._poly)
    def pos(self, an):
        return self._from_poly(an._poly)

    # -------------------
    # __float__
    # -------------------
    def float(self, an):
        return float(an._poly.eval(self.expr_float))

    
