from ifstype import *
from ifstype.rational import *
import timeit
from sympy import Rational as sy_R, sqrt
from sympy.polys.polytools import PurePoly
from sympy.abc import x
from ifstype.polynomial import Poly

from ifs_examples import *

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
    for name in finite_names:
        print(f"Running {name}")
        run_finite_type(examples.__dict__[name](),name)

def test_graph():
    g = Graph()
    for _ in range(6):
        g.add_vertex()
    g.add_edge(g.vertex(0), g.vertex(1))
    g.add_edge(g.vertex(1), g.vertex(2))
    g.add_edge(g.vertex(2), g.vertex(5))
    g.add_edge(g.vertex(5), g.vertex(1))
    g.add_edge(g.vertex(1), g.vertex(1))

    pos = sfdp_layout(g)
    comp, hist = label_components(g)
    graph_draw(g, pos, vertex_size=10, vertex_fill_color=comp, edge_pen_width=2,
               vcmap=matplotlib.cm.gist_heat_r, output="test.pdf")
    print(comp.a)
    #  print(comp)

if __name__ == "__main__":
    #  test_graph()
    run_finite_type(ifs6(),file_tag="test")
    #  run_infinite_type(ifs_osc(),stop=Rational(1,50),filename="infinite_non_overlap.pdf")
    #  get_alpha_density(ifs_osc(),count=20)

    #  run_infinite_type(ifs_osc(),stop=Rational(1,50),filename="infinite_non_overlap.pdf")
    #  get_alpha_density(ifs_osc(),count=20)



