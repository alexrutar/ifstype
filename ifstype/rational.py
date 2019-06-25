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
