from .numeric import Constants as C, Rational, AlgebraicNumberFactory
from .ifs import IFS,CtrFunc
from sympy import sqrt, pi
def ifs0():
    # finite type III with algebraic number
    expr = (sqrt(5)-1)/2
    anf = AlgebraicNumberFactory(expr)
    r = anf.alpha()
    return IFS.uniform_p(
            CtrFunc(r,Rational(0)),
            CtrFunc(r,Rational(1)-r))

def ifs1():
    # finite type III
    return IFS.uniform_p(
        CtrFunc(Rational(1,2),0),
        CtrFunc(Rational(1,4),Rational(1,4)),
        CtrFunc(Rational(1,2),Rational(1,2)))

def ifs2():
    # finite type IV
    return IFS.uniform_p(
        CtrFunc(Rational(1,3),0),
        CtrFunc(Rational(1,5),Rational(4,15)),
        CtrFunc(Rational(1,3),Rational(7,15)),
        CtrFunc(Rational(1,5),Rational(4,5)))
def ifs2m():
    # finite type IV
    return IFS.uniform_p(
        CtrFunc(Rational(1,3),0),
        CtrFunc(Rational(1,5),Rational(4,15)),
        CtrFunc(Rational(1,5),Rational(4,5)))

def ifs3():
    # not finite type
    return IFS.uniform_p(
        CtrFunc(Rational(1,3),0),
        CtrFunc(Rational(1,4),Rational(4,15)),
        CtrFunc(Rational(1,3),Rational(7,15)),
        CtrFunc(Rational(1,5),Rational(4,5)))

def ifs4():
    # finite type IV
    return IFS.uniform_p(
        CtrFunc(Rational(1,3),0),
        CtrFunc(Rational(1,4),Rational(1,4)),
        CtrFunc(Rational(1,3),Rational(1,2)),
        CtrFunc(Rational(1,4),Rational(3,4)))

def ifs5():
    # test this
    r = Rational(100,314)
    return IFS.uniform_p(
        CtrFunc(r,0),
        CtrFunc(r,r),
        CtrFunc(r,1))
