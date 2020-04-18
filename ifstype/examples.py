from .exact import Constants as C, NumberField, Poly, Fraction
from .ifs import IFS,AffineFunc, ifs_family


# ---------------------------------------------------------------------------
# Open Set Condition
# ---------------------------------------------------------------------------
@ifs_family
def osc_1(probs=None):
    return [AffineFunc(Fraction(1,3),Fraction(0)),
            AffineFunc(Fraction(1,3),Fraction(2,3))]

@ifs_family
def osc_2(probs=None):
    return [AffineFunc(Fraction(3,7),0),
            AffineFunc(Fraction(1,3),Fraction(2,3))]

@ifs_family
def osc_3(probs=None,a=None):
    # open set condition with overlap
    if a is None:
        a = Fraction(1,2)
    assert a <= Fraction(1,2)
    return [AffineFunc(a,0),
            AffineFunc(a*a,a-a*a),
            AffineFunc(a*a,1-a*a)]


# ---------------------------------------------------------------------------
# Equicontractive Finite Type
# ---------------------------------------------------------------------------
@ifs_family
def eft_1(probs=None):
    num_field = NumberField(Poly((-1,1,1)), 0.6180339887498948482045868)
    r = num_field.alpha()
    return [AffineFunc(r,Fraction(0)),
            AffineFunc(r,Fraction(1)-r)]


@ifs_family
def eft_k(probs=None):
    r = Fraction(1,3)
    return [AffineFunc(r,0),
            AffineFunc(r,1*r),
            AffineFunc(r,2*r),
            AffineFunc(r,Fraction(1,9))]

@ifs_family
def eft_2(probs=None,r=Fraction(1,5)):
    assert 0 <  r <= Fraction(1,2)
    return [AffineFunc(r,0),
            AffineFunc(r,r),
            AffineFunc(r,1)]

@ifs_family
def eft_3(probs=None):
    num_field = NumberField(Poly((-1,1,1)), 0.6180339887498948482045868)
    r = num_field.alpha()
    return [AffineFunc(r,0),
            AffineFunc(r,Fraction(1,2)-r*Fraction(1,2)),
            AffineFunc(r,Fraction(1)-r)]

@ifs_family
def eft_4(probs=None):
    num_field = NumberField(Poly((-1,1,0,1)),0.6823278038280193273694837)
    r = num_field.alpha()
    return [AffineFunc(r,0),
            AffineFunc(r,1-r)]

@ifs_family
def eft_5(probs=None):
    num_field = NumberField(Poly((-1,-1,-1,1)), 1.83928675521416113255185)
    r = 1/num_field.alpha()
    return [AffineFunc(r,0),
            AffineFunc(r,r*r),
            AffineFunc(r,r)]

@ifs_family
def eft_6(probs=None):
    r = Fraction(1,3)
    return [AffineFunc(r,0),
            AffineFunc(r,Fraction(2,87)),
            AffineFunc(r,Fraction(2,3))]

# ---------------------------------------------------------------------------
# Finite Type
# ---------------------------------------------------------------------------
@ifs_family
def ft_1(probs=None):
    return [AffineFunc(Fraction(1,2),0),
            AffineFunc(Fraction(1,4),Fraction(1,4)),
            AffineFunc(Fraction(1,2),Fraction(1,2))]

@ifs_family
def ft_2(probs=None):
    num_field = NumberField(Poly((-1,1,1)), 0.6180339887498948482045868)
    r = num_field.alpha()
    return [AffineFunc(r,0),
            AffineFunc(-r,1)]
@ifs_family
def ft_3(probs=None,a=None):
    if a is None:
        a = Fraction(1,3)
    # open set condition with overlap
    return [AffineFunc(a,0),
            AffineFunc(a*a,a-a*a),
            AffineFunc(a,1-a)]
# ---------------------------------------------------------------------------
# Finite Neighbour Condition
# ---------------------------------------------------------------------------
@ifs_family
def fnc_1(probs=None):
    "Example with full essential class."
    return [AffineFunc(Fraction(1,3),0),
            AffineFunc(Fraction(1,5),Fraction(4,15)),
            AffineFunc(Fraction(1,3),Fraction(7,15)),
            AffineFunc(Fraction(1,5),Fraction(4,5))]


@ifs_family
def fnc_2(probs=None,a=Fraction(1,4),b=Fraction(1,3)):
    "Example with no return to 0"
    assert a+b-a*b <= Fraction(1,2) and 0 < a and 0 < b
    return [AffineFunc(a,0),
            AffineFunc(b,a-a*b),
            AffineFunc(a,1-a-b+a*b),
            AffineFunc(b,1-b)]

@ifs_family
def fnc_3(probs=None,a=None):
    "Example with no return to 0 and single element loop class"
    # can also remove last function and is non-trivial IFS subset
    # example where reduced neighbours are not the same as neighbours
    # also has, at first pass, a neighbour set with 0 children
    num_field = NumberField(Poly((-1,3,2)),0.28077640640441513745535246399) # upper bound on arbitrary params
    r = num_field.alpha()
    if a is None:
        a=r
    assert a <= r
    b = a/(1+a)
    return [AffineFunc(a,0),
            AffineFunc(b,a*b),
            AffineFunc(a,1-Fraction(2)*b),
            AffineFunc(b,1-b)]

@ifs_family
def fnc_4(probs=None,a=None):
    """Example with loop class at 0 and remainder essential class.
    When a is maximal, this IFS is in fact finite type and has 3 loop class and essential class size 4"""
    num_field = NumberField(Poly((-1,1,1)), 0.6180339887498948482045868)
    r = num_field.alpha()
    if a is None:
        a=r
    assert 0<a <= r
    b = a/(1+a)
    return [AffineFunc(a,0),
            AffineFunc(b,a*b),
            AffineFunc(b,1-b)]


@ifs_family
def fnc_5(probs=None,a=None,b=None):
    if a is None:
        a = Fraction(1,3)
    if b is None:
        b = Fraction(1,4)
    assert a+2*b-a*b <= 1
    return [AffineFunc(a,0),
            AffineFunc(b,a*(1-b)),
            AffineFunc(b,1-b)]

# some new experimental examples
@ifs_family
def fnc_6(probs=None,a=None,b=None,n=None):
    if a is None:
        a = Fraction(1,4)
    if b is None:
        b = Fraction(1,3)
    if n is None:
        n = 3
    assert a+b-a*(b**n)+b <= 1
    return [AffineFunc(a,0),
            AffineFunc(b**n,a-a*(b**n)),
            AffineFunc(b,1-b)]

@ifs_family
def fnc_7(probs=None):
    a = Fraction(1,4)
    b = Fraction(1,3)
    n = 50
    return [AffineFunc(a,0),
            AffineFunc(b**n,a-a*(b**n)),
            AffineFunc(a,1-a-b+a*b),
            AffineFunc(b,1-b)]

@ifs_family
def fnc_8(probs=None):
    a = Fraction(1,3)
    b = Fraction(1,4)
    return [AffineFunc(a,0),
            AffineFunc(b,a-a*b),
            AffineFunc(b,2*b),
            AffineFunc(b,3*b)]

# ---------------------------------------------------------------------------
# Not Finite Neighbour Condition
# ---------------------------------------------------------------------------
@ifs_family
def inf_1(probs=None):
    return [AffineFunc(Fraction(3,7),0),
            AffineFunc(Fraction(2,3),Fraction(1,3))]

@ifs_family
def inf_2(probs=None):
    return [AffineFunc(Fraction(1,2),0),
            AffineFunc(Fraction(1,5),Fraction(3,10)),
            AffineFunc(Fraction(1,5),Fraction(4,5))]

@ifs_family
def inf_3(probs=None):
    # no return to 0 and single element loop class
    return [AffineFunc(Fraction(1,3),0),
            AffineFunc(Fraction(1,4),Fraction(1,6)),
            AffineFunc(Fraction(1,3),Fraction(3,5)),
            AffineFunc(Fraction(1,4),Fraction(4,5))]

@ifs_family
def inf_4(probs=None):
    return [AffineFunc(Fraction(1,7),0),
            AffineFunc(Fraction(1,7),Fraction(2,21)),
            AffineFunc(Fraction(1,7),Fraction(4,21)),
            AffineFunc(Fraction(1,7),Fraction(6,7))]

@ifs_family
def test(probs=None):
    return [AffineFunc(Fraction(1,3),0),
            AffineFunc(Fraction(1,4),Fraction(4,5)),
            AffineFunc(-Fraction(1,2),1)]

@ifs_family
def test2(probs=None):
    return [AffineFunc(Fraction(1,3),0),
            AffineFunc(Fraction(1,4),Fraction(1,5)),
            AffineFunc(Fraction(3,10),Fraction(1,5)+Fraction(1,4)),
            AffineFunc(Fraction(1,4),Fraction(3,4))]

# ---------------------------------------------------------------------------
# Main Examples
# ---------------------------------------------------------------------------
@ifs_family
def example_osc(probs=None,a=None):
    # open set condition with overlap
    if a is None:
        a = Fraction(1,2)
    assert a <= Fraction(1,2)
    return [AffineFunc(a,0),
            AffineFunc(a*a,a-a*a),
            AffineFunc(a*a,1-a*a)]

@ifs_family
def example_eft(probs=None):
    return [AffineFunc(Fraction(1,3),0),
            AffineFunc(Fraction(1,3),Fraction(1,9)),
            AffineFunc(Fraction(1,3),Fraction(1,3)),
            AffineFunc(Fraction(1,3),Fraction(2,3))]

@ifs_family
def example_fnc(probs=None):
    return [AffineFunc(Fraction(1,4),0),
            AffineFunc(Fraction(1,3),Fraction(1,6)),
            AffineFunc(Fraction(1,4),Fraction(1,2)),
            AffineFunc(Fraction(1,4),Fraction(3,4))]

