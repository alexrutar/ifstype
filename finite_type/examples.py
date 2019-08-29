from sympy import sqrt, Rational as sy_R
import inspect

from .numerics.rational import Constants as C, Rational
from .numerics.algebraic import AlgebraicNumberFactory
from .numerics.polynomial import Poly
from .ifs import IFS,CtrFunc

def ifs_family(fn):
    """Decorator to construct families of iterated function systems parametrized by some set of functions.
    When decorated with @ifs_family, call functions by first specifying probs if necessary, followed by keyword arguments.
    """
    def wrapper(**kwargs):
        fn_kwargs = {k:v.default for k,v in inspect.signature(fn).parameters.items()}
        func_params = {**fn_kwargs,**kwargs}
        try:
            probs = func_params['probs']
        except KeyError:
            raise KeyError("Function decorated by 'ifs_family' has no default keyword 'probs'.")
        if probs is None:
            return IFS.uniform_p(*fn(**func_params))
        else:
            return IFS(fn(**func_params),probs)
    return wrapper


# ---------------------------------------------------------------------------
# Open Set Condition
# ---------------------------------------------------------------------------
@ifs_family
def osc_1(probs=None):
    return [CtrFunc(Rational(1,3),Rational(0)),
            CtrFunc(Rational(1,3),Rational(2,3))]

@ifs_family
def osc_2(probs=None):
    return [CtrFunc(Rational(3,7),0),
            CtrFunc(Rational(1,3),Rational(2,3))]

@ifs_family
def osc_3(probs=None):
    # open set condition with overlap
    return [CtrFunc(Rational(1,2),0),
            CtrFunc(Rational(1,4),Rational(1,4)),
            CtrFunc(Rational(1,4),Rational(1,4))]


# ---------------------------------------------------------------------------
# Equicontractive Finite Type
# ---------------------------------------------------------------------------
@ifs_family
def eft_1(probs=None):
    # finite type III with algebraic number
    expr = (sqrt(5)-1)/2
    anf = AlgebraicNumberFactory.from_sympy(expr)
    r = anf.alpha()
    return [CtrFunc(r,Rational(0)),
            CtrFunc(r,Rational(1)-r)]

@ifs_family
def eft_2(probs=None,r=Rational(1,5)):
    assert 0 <  r <= Rational(1,2)
    return [CtrFunc(r,0),
            CtrFunc(r,r),
            CtrFunc(r,1)]

@ifs_family
def eft_3(probs=None):
    expr = sy_R(1,2)*(sqrt(5)-sy_R(1))
    anf = AlgebraicNumberFactory.from_sympy(expr)
    r = anf.alpha()
    return [CtrFunc(r,0),
            CtrFunc(r,Rational(1,2)-r*Rational(1,2)),
            CtrFunc(r,Rational(1)-r)]

@ifs_family
def eft_4(probs=None):
    anf = AlgebraicNumberFactory(Poly((-1,1,0,1)),0.6823278038280193273694837)
    r = anf.alpha()
    return [CtrFunc(r,0),
            CtrFunc(r,1-r)]

@ifs_family
def eft_5(probs=None):
    anf = AlgebraicNumberFactory(Poly((-1,-1,-1,1)), 1.83928675521416113255185)
    r = 1/anf.alpha()
    return [CtrFunc(r,0),
            CtrFunc(r,r*r),
            CtrFunc(r,r)]

# ---------------------------------------------------------------------------
# Finite Type
# ---------------------------------------------------------------------------
@ifs_family
def ft_1(probs=None):
    return [CtrFunc(Rational(1,2),0),
            CtrFunc(Rational(1,4),Rational(1,4)),
            CtrFunc(Rational(1,2),Rational(1,2))]

@ifs_family
def ft_2(probs=None):
    anf = AlgebraicNumberFactory(Poly((-1,1,1)), 0.6180339887498948482045868)
    r = anf.alpha()
    return [CtrFunc(r,0),
            CtrFunc(-r,1)]
# ---------------------------------------------------------------------------
# Weak Finite Type
# ---------------------------------------------------------------------------
@ifs_family
def wft_1(probs=None):
    "Example with full essential class."
    return [CtrFunc(Rational(1,3),0),
            CtrFunc(Rational(1,5),Rational(4,15)),
            CtrFunc(Rational(1,3),Rational(7,15)),
            CtrFunc(Rational(1,5),Rational(4,5))]


@ifs_family
def wft_2(probs=None,a=Rational(1,4),b=Rational(1,3)):
    "Example with no return to 0"
    assert a+b-a*b <= Rational(1,2) and 0 < a and 0 < b
    return [CtrFunc(a,0),
            CtrFunc(b,a-a*b),
            CtrFunc(a,1-a-b+a*b),
            CtrFunc(b,1-b)]

@ifs_family
def wft_3(probs=None):
    "Example with no return to 0 and single element loop class"
    return [CtrFunc(Rational(1,4),0),
            CtrFunc(Rational(1,5),Rational(1,20)),
            CtrFunc(Rational(1,4),Rational(3,5)),
            CtrFunc(Rational(1,5),Rational(4,5))]

@ifs_family
def wft_4(probs=None,a=None):
    "Example with 3 loop classes and small essential class (size 4)"
    anf = AlgebraicNumberFactory(Poly((-1,1,1)), 0.6180339887498948482045868)
    r = anf.alpha()
    if a is None:
        a=r
    assert 0<a <= r
    b = a/(1+a)
    return [CtrFunc(a,0),
            CtrFunc(b,a*b),
            CtrFunc(b,1-b)]


# ---------------------------------------------------------------------------
# Not Weak Finite Type
# ---------------------------------------------------------------------------
@ifs_family
def inf_1(probs=None):
    return [CtrFunc(Rational(3,7),0),
            CtrFunc(Rational(2,3),Rational(1,3))]

@ifs_family
def inf_2(probs=None):
    return [CtrFunc(Rational(1,2),0),
            CtrFunc(Rational(1,5),Rational(3,10)),
            CtrFunc(Rational(1,5),Rational(4,5))]

@ifs_family
def inf_3(probs=None):
    # no return to 0 and single element loop class
    return [CtrFunc(Rational(1,3),0),
            CtrFunc(Rational(1,4),Rational(1,6)),
            CtrFunc(Rational(1,3),Rational(3,5)),
            CtrFunc(Rational(1,4),Rational(4,5))]

