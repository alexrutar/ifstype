from ifstype import *
from ifstype.rational import *
import timeit
from sympy import Rational as sy_R, sqrt
from sympy.polys.polytools import PurePoly
from sympy.abc import x
from ifstype.polynomial import Poly

import ifs_examples

from graph_tool.all import *
import matplotlib

def test_poly_gcd():
    pol = PurePoly(x+sy_R(1,3))-PurePoly(x)
    pol2 = PurePoly(x**11)

    in_pol = Poly.from_PurePoly(pol)
    in_pol2 = Poly.from_PurePoly(pol2)

    for p1, p2 in zip(pol.gcdex(pol2), in_pol.gcdex(in_pol2)):
        print(f"{p1} : {p2}")

def test_poly_invert():
    pol = PurePoly(x+sy_R(1,3))-PurePoly(x)
    pol2 = PurePoly(x**11)

    in_pol = Poly.from_PurePoly(pol)
    in_pol2 = Poly.from_PurePoly(pol2)

    print(pol.invert(pol2))
    print(in_pol.invert(in_pol2))

def test_algebraics():
    expr = sy_R(1,3)*(1 + (19-3*sqrt(33))**sy_R(1,3) + (19+3*sqrt(33))**sy_R(1,3))
    fl = float(expr)
    anf = AlgebraicNumberFactory.from_sympy(expr)
    print(float(anf.alpha()/(Rational(1)+anf.alpha())))
    print(fl/(1+fl))

def run_all_finite():
    for name in ifs_examples.finite_names:
        print(f"Running {name}")
        run_finite_type(ifs_examples.__dict__[name](),"out/"+name)

if __name__ == "__main__":
    #  run_all_finite()
    #  test_graph()
    #  run_finite_type(ifs_examples.test_3(),file_tag="test")
    ifs_examples.test_1()
    ifs_examples.test_2()
    ifs_examples.test_3()
    #  run_infinite_type(ifs_osc(),stop=Rational(1,50),filename="infinite_non_overlap.pdf")
    #  get_alpha_density(ifs_osc(),count=20)

    #  run_infinite_type(ifs_osc(),stop=Rational(1,50),filename="infinite_non_overlap.pdf")
    #  get_alpha_density(ifs_osc(),count=20)



