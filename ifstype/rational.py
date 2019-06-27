from quicktions import Fraction

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

Constants = _Const()

class _Inf:
    pass

class _PosInf(_Inf):
    def __eq__(self,other):
        if isinstance(other,PosInf):
            return True
        return False
    def __neq__(self,other):
        if isinstance(other,PosInf):
            return False
        return True

    def __lt__(self,other):
        return False
    def __le__(self,other):
        return self == other or self < other
    def __gt__(self,other):
        return True
    def __ge__(self,other):
        return self == other or self > other

    def __add__(self,other):
        if isinstance(other,_NegInf):
            raise ValueError("PosInf + NegInf is undefined")
        else:
            return PosInf

    def __radd__(self,other):
        if isinstance(other,_NegInf):
            raise ValueError("PosInf + NegInf is undefined")
        else:
            return PosInf

    def __sub__(self,other):
        if isinstance(other,_PosInf):
            raise ValueError("PosInf - PosInf is undefined")
        else:
            return PosInf

    def __rsub__(self,other):
        if isinstance(other,_PosInf):
            raise ValueError("PosInf - PosInf is undefined")
        else:
            return NegInf

    def __mul__(self,other):
        if other == 0:
            raise ValueError("PosInf * 0 is undefined")
        elif other > 0:
            return PosInf
        else:
            return NegInf

    def __rmul__(self,other):
        if other == 0:
            raise ValueError("PosInf * 0 is undefined")
        elif other > 0:
            return PosInf
        else:
            return NegInf

    def __div__(self,other):
        if other == 0:
            raise ZeroDivisionError("PosInf / 0 is undefined")
        elif isinstance(other,_Inf):
            raise ValueError("Inf / Inf is undefined")
        elif other > 0:
            return PosInf
        else:
            return NegInf

    def __rdiv__(self,other):
        raise ValueError("Inf / Inf is undefined")

    def __neg__(self):
        return NegInf

    def __pos__(self):
        return PosInf

    def __abs__(self):
        return PosInf

class _NegInf(_Inf):
    def __eq__(self,other):
        if isinstance(other,NegInf):
            return True
        return False
    def __neq__(self,other):
        if isinstance(other,NegInf):
            return False
        return True

    def __lt__(self,other):
        return True
    def __le__(self,other):
        return self == other or self < other
    def __gt__(self,other):
        return False
    def __ge__(self,other):
        return self == other or self > other

    def __add__(self,other):
        if isinstance(other,_PosInf):
            raise ValueError("PosInf + NegInf is undefined")
        else:
            return NegInf

    def __radd__(self,other):
        if isinstance(other,_PosInf):
            raise ValueError("PosInf + NegInf is undefined")
        else:
            return NegInf

    def __sub__(self,other):
        if isinstance(other,_NegInf):
            raise ValueError("NegInf - NegInf is undefined")
        else:
            return NegInf

    def __rsub__(self,other):
        if isinstance(other,_NegInf):
            raise ValueError("NegInf - NegInf is undefined")
        else:
            return PosInf

    def __mul__(self,other):
        if other == 0:
            raise ValueError("NegInf * 0 is undefined")
        elif other > 0:
            return NegInf
        else:
            return PosInf

    def __rmul__(self,other):
        if other == 0:
            raise ValueError("PosInf * 0 is undefined")
        elif other > 0:
            return PosInf
        else:
            return NegInf

    def __div__(self,other):
        if other == 0:
            raise ZeroDivisionError("PosInf / 0 is undefined")
        elif isinstance(other,_Inf):
            raise ValueError("Inf / Inf is undefined")
        elif other > 0:
            return NegInf
        else:
            return PosInf

    def __rdiv__(self,other):
        raise ValueError("Inf / Inf is undefined")

    def __neg__(self):
        return PosInf

    def __pos__(self):
        return NegInf

    def __abs__(self):
        return PosInf

PosInf = _PosInf()
NetInf = _NegInf()
