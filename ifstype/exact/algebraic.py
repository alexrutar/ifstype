import math
import itertools
import functools
import numbers

from .polynomial import Poly
from .rational import Fraction, Constants as C

# parse the args to convert a Fraction argument to AlgebraicNumber, to be used with NumberField
def check_other(f):
    "Decorator to verify that algebraic numbers originate from equivalent factories, and convertes Fractions to correct format."
    def valid_num_field(self,other):
        if isinstance(other,numbers.Rational):
            return f(self, self.num_field._from_poly(Poly.cnst(other)))
        elif isinstance(other,AlgebraicNumber) and self.num_field == other.num_field:
            return f(self, other)
        else:
            raise ValueError(f"Inequivalent algebraic numbers! Other is a {type(other)}.")
    return valid_num_field

class AlgebraicNumber(numbers.Real):
    def __init__(self, num_field, poly):
        self.num_field = num_field # pointer to the factory
        self._poly = poly

    def __repr__(self):
        return f"AlgebraicNumber(expr={self._poly}, x={self.num_field.expr_float})"

    def __str__(self):
        return f"{self._poly.with_symbol(self.num_field.symbol,space=False)}"

    def as_str(self, symbol=None, term=False, space=False):
        "String method, also compatible with latex."
        if symbol is None:
            symbol = self.num_field.symbol
        base_str = self._poly.with_symbol(symbol,space=space)
        if term:
            num_terms = len([x for x in self._poly if x != C.n_0])
            if num_terms > 1:
                base_str = f"({base_str})"
        return base_str

    def as_latex(self,**kwargs):
        return self.as_str(**kwargs)

    def __hash__(self):
        return hash(self._poly)

    # -----------------------------------------
    # comparisons
    # -----------------------------------------
    @check_other
    def __lt__(self, other):
        return self.num_field.lt(self,other)

    @check_other
    def __le__(self, other):
        return self.num_field.le(self,other)

    @check_other
    def __gt__(self, other):
        return self.num_field.gt(self,other)

    @check_other
    def __ge__(self, other):
        return self.num_field.ge(self,other)

    @check_other
    def __eq__(self,other):
        return self.num_field.eq(self,other)

    @check_other
    def __neq__(self,other):
        return self.num_field.neq(self,other)
    # -----------------------------------------
    # numeric
    # -----------------------------------------

    @check_other
    def __truediv__(self, other):
        return self.num_field.mul(self,self.num_field.inv(other))

    @check_other
    def __rtruediv__(self, other):
        return self.num_field.mul(self.num_field.inv(self),other)

    @check_other
    def __add__(self,other):
        return self.num_field.add(self,other)

    @check_other
    def __radd__(self,other):
        return self.num_field.add(self,other)

    @check_other
    def __sub__(self,other):
        return self.num_field.sub(self,other)

    @check_other
    def __rsub__(self,other):
        return self.num_field.sub(other,self)

    @check_other
    def __mul__(self,other):
        return self.num_field.mul(self,other)

    @check_other
    def __rmul__(self,other):
        return self.num_field.mul(self,other)
        
    def __pow__(self, it):
        return self.num_field.pow(self,it)

    # ---------------------------------
    # unary operations
    # ---------------------------------
    def __abs__(self):
        return self.num_field.abs(self)
    def __neg__(self):
        return self.num_field.neg(self)
    def __pos__(self):
        return self.num_field.pos(self)

    # ---------------------------------
    # support for other numeric operations
    # ---------------------------------
    def __float__(self):
        return self.num_field.float(self)


    def __floor__(self):
        return float(self).__floor__()
    def __ceil__(self):
        return float(self).__ceil__()

    def __floordiv__(self,other):
        return float(self).__floordiv__(other)
    def __rfloordiv__(self,other):
        return float(self).__rfloordiv__(other)

    def __mod__(self,other):
        return float(self).__mod__(other)
    def __rmod__(self,other):
        return float(self).__rmod__(other)

    def __round__(self, ndigits=None):
        return float(self).__round__(self, ndigits=ndigits)
    def __trunc__(self):
        return float(self).__trunc__()

    def __rpow__(self, base):
        return float(self).__rpow__(base)

# an algebraic number is a rational
AlgebraicNumber.register(numbers.Rational)

class NumberField:
    """A class governing the creation of algebraic numbers.
    Compute the number field Q(alpha) where alpha ~ expr_float has minimal polynomial `minpoly` (a exact.Poly instance).
    Two factories are considered equivalent if they come from the same expression.
    Represents algebraic numbers with respect to their basis over Q.
    Internally, for comparisons, the float form is used in the expression. Be careful, if you need a lot of precision, this could be dangerous!"""
    # TODO: have arbitrary precision float / detect when floats may not be precise enough.
    # for most applications, this should be fine, at least in finite type IFS since all the comparisons will be with intervals about the same length
    def __init__(self, minpoly, expr_float, symbol='a'):
        self.minpoly = minpoly
        self.expr_float = expr_float
        self.deg = self.minpoly.deg
        self.symbol = symbol

    @classmethod
    def from_sympy(cls, sy_expr):
        from sympy.polys.domains import QQ
        from sympy.polys.numberfields import minpoly

        minpoly = Poly.from_PurePoly(minpoly(sy_expr,polys=True,domain=QQ))
        expr_float = float(sy_expr)
        return cls(minpoly,expr_float)

    #  def an_string_rep(self, an):
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

    
