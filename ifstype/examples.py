from .exact import Constants as C, Rational, NumberField, Poly
from .ifs import IFS,CtrFunc, ifs_family


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
    num_field = NumberField(Poly((-1,1,1)), 0.6180339887498948482045868)
    r = num_field.alpha()
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
    num_field = NumberField(Poly((-1,1,1)), 0.6180339887498948482045868)
    r = num_field.alpha()
    return [CtrFunc(r,0),
            CtrFunc(r,Rational(1,2)-r*Rational(1,2)),
            CtrFunc(r,Rational(1)-r)]

@ifs_family
def eft_4(probs=None):
    num_field = NumberField(Poly((-1,1,0,1)),0.6823278038280193273694837)
    r = num_field.alpha()
    return [CtrFunc(r,0),
            CtrFunc(r,1-r)]

@ifs_family
def eft_5(probs=None):
    num_field = NumberField(Poly((-1,-1,-1,1)), 1.83928675521416113255185)
    r = 1/num_field.alpha()
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
    num_field = NumberField(Poly((-1,1,1)), 0.6180339887498948482045868)
    r = num_field.alpha()
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
def wft_3(probs=None,a=None):
    "Example with no return to 0 and single element loop class"
    num_field = NumberField(Poly((-1,3,2)),0.28077640640441513745535246399) # upper bound on arbitrary params
    r = num_field.alpha()
    if a is None:
        a=r
    assert a <= r
    b = a/(1+a)
    return [CtrFunc(a,0),
            CtrFunc(b,a*b),
            CtrFunc(a,1-Rational(2)*b),
            CtrFunc(b,1-b)]

@ifs_family
def wft_4(probs=None,a=None):
    """Example with loop class at 0 and remainder essential class.
    When a is maximal, this IFS is in fact finite type and has 3 loop class and essential class size 4"""
    num_field = NumberField(Poly((-1,1,1)), 0.6180339887498948482045868)
    r = num_field.alpha()
    if a is None:
        a=r
    assert 0<a <= r
    b = a/(1+a)
    return [CtrFunc(a,0),
            CtrFunc(b,a*b),
            CtrFunc(b,1-b)]

@ifs_family
def wft_5(probs=None,a=None):
    """Example with loop class at 0 and remainder essential class.
    When a is maximal, this IFS is in fact finite type and has 3 loop class and essential class size 4"""
    num_field = NumberField(Poly((-1,1,1)), 0.6180339887498948482045868)
    r = num_field.alpha()
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

