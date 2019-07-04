from ifstype.rational import Constants as C, Rational
from ifstype.algebraic import AlgebraicNumberFactory
from ifstype.polynomial import Poly
from ifstype.ifs import IFS,CtrFunc
from sympy import sqrt, Rational as sy_R

finite_names = ['ifs0','ifs1','ifs2','ifs3','ifs4','ifs5','ifs6','ifs7']

def ifs0():
    # finite type III with algebraic number
    expr = (sqrt(5)-1)/2
    anf = AlgebraicNumberFactory.from_sympy(expr)
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

def ifs3():
    # finite type IV
    return IFS.uniform_p(
        CtrFunc(Rational(1,3),0),
        CtrFunc(Rational(1,4),Rational(1,4)),
        CtrFunc(Rational(1,3),Rational(1,2)),
        CtrFunc(Rational(1,4),Rational(3,4)))

def ifs4():
    # probably equicontractive finite type
    r = Rational(111,371)
    return IFS.uniform_p(
        CtrFunc(r,0),
        CtrFunc(r,r),
        CtrFunc(r,1))

def ifs5(rho=Rational(1,5),r=Rational(3,7)):
    # very similar to ifs2
    assert rho+2*r-rho*r <= 1
    return IFS.uniform_p(
            CtrFunc(rho,0),
            CtrFunc(r,rho*(1-r)),
            CtrFunc(r,1-r))

def ifs6():
    expr = sy_R(1,2)*(sqrt(5)-sy_R(1))
    anf = AlgebraicNumberFactory.from_sympy(expr)
    r = anf.alpha()
    return IFS.uniform_p(
            CtrFunc(r,0),
            CtrFunc(r,Rational(1,2)-r*Rational(1,2)),
            CtrFunc(r,Rational(1)-r))

def ifs7():
    anf = AlgebraicNumberFactory(Poly((-1,1,0,1)),0.6823278038280193273694837)
    r = anf.alpha()
    return IFS.uniform_p(
            CtrFunc(r,0),
            CtrFunc(r,1-r))

def test_0():
    n = 3
    z_list = [-1,3,4,7]
    return IFS.uniform_p(*[CtrFunc(Rational(1,n),z) for z in z_list])

def test_1():
    # probably equicontractive finite type
    r = Rational(111,371)
    return IFS.uniform_p(
        CtrFunc(r,0),
        CtrFunc(r,r),
        CtrFunc(r,1))

def test_2(rho=Rational(1,5),r=Rational(3,7)):
    # very similar to ifs2
    assert rho+2*r-rho*r <= 1
    return IFS.uniform_p(
            CtrFunc(rho,0),
            CtrFunc(r,rho*(1-r)),
            CtrFunc(r,1-r))

def test_3():
    expr = sy_R(1,3)*(1 + (19-3*sqrt(33))**sy_R(1,3) + (19+3*sqrt(33))**sy_R(1,3))
    anf = AlgebraicNumberFactory.from_sympy(expr)
    r = 1/anf.alpha()
    return IFS.uniform_p(
            CtrFunc(r,0),
            CtrFunc(r,r**2),
            CtrFunc(r,r))

def test_4(rho=Rational(1,5),r=Rational(3,7)):
    # very similar to ifs2
    assert rho+2*r-rho*r <= 1
    return IFS.uniform_p(
            CtrFunc(rho,0),
            CtrFunc(r,rho*(1-r)),
            CtrFunc(r,1-r))


# infinite type
def inf_ifs0():
    return IFS.uniform_p(
            CtrFunc(Rational(3,7),0),
            CtrFunc(Rational(2,3),Rational(1,3)))

def ifs_osc():
    # basic equicontractive ifs with no overlap
    return IFS.uniform_p(
            CtrFunc(Rational(3,7),0),
            CtrFunc(Rational(1,3),Rational(2,3)))


def inf_ifs3():
    # not finite type
    return IFS.uniform_p(
        CtrFunc(Rational(1,3),0),
        CtrFunc(Rational(1,4),Rational(4,15)),
        CtrFunc(Rational(1,3),Rational(7,15)),
        CtrFunc(Rational(1,5),Rational(4,5)))
