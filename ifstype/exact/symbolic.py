import itertools
import numpy as np
from operator import mul
from functools import reduce
from collections import defaultdict
from quicktions import Fraction

from .algebraic import AlgebraicNumber

# TODO: implement __radd__ etc. and do lookup for rational numbers, to allow scalar multiplication and adding rationals
def check_eq(fn):
    def wrapper(self, other):
        if isinstance(other,AlgebraicNumber):
            return fn(self, self.ring.from_algebraic(other))
        elif isinstance(other, SymbolicElement):
            if not self.ring is other.ring:
                raise ValueError("Symbolic elements originate from distinct SymbolicRing instances")
            return fn(self,other)
        else:
            raise ValueError("Operations between types not defined!")

    return wrapper

class _PowerIdx(tuple):
    def __new__(cls, idx):
        self = super().__new__(cls,idx)
        return self

    def __add__(self, other):
        return _PowerIdx(s+o for s,o in zip(self,other))

class SymbolicRing:
    """Governing ring of symbolic elements on a fixed set of symbols."""
    def __init__(self, symbols, evals=None):
        self.symb_dct = {symbol:idx for idx,symbol in enumerate(symbols)}
        self.evals = evals

    def get_symbols(self):
        return tuple(self.term(s) for s in self.symb_dct.keys())

    def term(self, symbs,coef=1):
        new_dct = {_PowerIdx(1 if p in symbs else 0 for p in self.symb_dct.keys()):coef}
        return SymbolicElement(new_dct, self)

    def from_algebraic(self, rt):
        new_dct = {_PowerIdx(0 for _ in range(len(self.symb_dct.keys()))):rt}
        return SymbolicElement(new_dct, self)
    
    def from_dct(self, dct):
        return SymbolicElement(dct, self)

    def default_getter(self):
        return 0

    def from_dct(self, new_dct):
        return SymbolicElement({k:v for k,v in new_dct.items() if v != 0}, self)

    def is_zero(self, se):
        return len(se.coef_dct.keys()) == 0

    def eq(self, se1, se2):
        return self.is_zero(se1-se2)

    def add(self, se1, se2):
        "add two symbolic elements together"
        new_dct = defaultdict(self.default_getter,se1.coef_dct.copy())
        for k,v in se2.coef_dct.items():
            new_dct[k] += v
        return self.from_dct(new_dct)

    def sub(self, se1, se2):
        new_dct = defaultdict(self.default_getter,se1.coef_dct.copy())
        for k,v in se2.coef_dct.items():
            new_dct[k] -= v
        return self.from_dct(new_dct)

    def mul(self, se1, se2):
        new_dct = defaultdict(self.default_getter)
        for (id1,v1),(id2,v2) in itertools.product(se1.coef_dct.items(), se2.coef_dct.items()):
            new_dct[id1+id2] += v1*v2
        return self.from_dct(new_dct)

    def pow(self, se, it):
        if it < 0 or not isinstance(it,int):
            raise NotImplementedError("Power must be integer >= 0")
        elif it == 0:
            return self.from_algebraic(1)
        elif it == 1:
            return self.from_dct(se.coef_dct)
        elif it == 2:
            return self.mul(se,se)
        # use divide and conquer algorithm for larger powers
        else:
            return self.mul(self.pow(se,-(-it // 2)),self.pow(se,(it // 2)))

    def sig_to_str(self, sig):
        reverse = {v:idx for idx,v in self.symb_dct.items()}
        def check_val(v):
            if v == 1:
                return ""
            else:
                return f"^{v}"

        return "*".join(f"{reverse[idx]}{check_val(v)}" for idx,v in enumerate(sig) if v != 0)
    
    def set_eval(self, val_dct):
        assert all(symb in val_dct.keys() for symb in self.symb_dct.keys()), "Valuation dict must specify values for all variables."
        self.evals = {self.symb_dct[k]:v for k,v in val_dct.items()}

    def _eval_term(self, coef,sig):
        if self.evals is not None:
            return coef * reduce(mul,(v**sig[k] for k,v in self.evals.items()))
        else:
            raise ValueError("Must set eval before float conversion!")


class SymbolicElement:
    def __init__(self, coef_dct, ring):
        self.ring = ring
        self.coef_dct = coef_dct

    def __hash__(self):
        return hash(str(self))

    def __float__(self):
        return float(self.eval())

    def eval(self):
        "An eval_dct is a dictionary of {symb:val} to substitute into the string"
        return sum(self.ring._eval_term(v, sig) for sig,v in self.coef_dct.items())

    @check_eq
    def __eq__(self, other):
        return self.ring.eq(self,other)
    @check_eq
    def __neq__(self, other):
        return not self.ring.sq(self,other)

    @check_eq
    def __add__(self, other):
        return self.ring.add(self,other)
    @check_eq
    def __radd__(self,other):
        return self.ring.add(self, other)

    @check_eq
    def __mul__(self, other):
        return self.ring.mul(self,other)
    @check_eq
    def __rmul__(self, other):
        return self.ring.mul(self,other)

    @check_eq
    def __sub__(self, other):
        return self.ring.sub(self,other)

    @check_eq
    def __rsub__(self,other):
        return self.ring.sub(other,self)

    def __pow__(self, it):
        return self.ring.pow(self,it)

    def __str__(self):
        terms = "".join((" + " if val >= 0 else " - ") + (f"{abs(val)}*{self.ring.sig_to_str(sig)}" if val != 1 else self.ring.sig_to_str(sig)) for sig,val in self.coef_dct.items() if val != 0)
        if len(terms) == 0:
            return "0"
        elif terms[1] == "+":
            return terms[3:]
        else:
            return terms[1:]

class SymbolicMatrix:
    """TransitionMatrix class to represent the transition matrix associated to an edge in the transition graph.
    """
    def __init__(self, double_list):
        self.matrix = tuple(tuple(sublist) for sublist in double_list)
        self.xdim = len(self.matrix)
        self.ydim = len(self.matrix[0])
        if not all(len(self.matrix[i]) == self.ydim for i in range(self.xdim)):
            raise ValueError("Invalid argument dimensions")

    @classmethod
    def identity(cls,n):
        dbl_lst = [[0 for _ in range(n)] for _ in range(n)]
        for i in range(n):
            dbl_list[i][i] = 1
        return cls(dbl_list)

    def __array__(self):
        return np.array(self.matrix, dtype=float)

    def __hash__(self):
        return hash(self.matrix.tobytes())

    def __repr__(self):
        return repr(self.matrix)

    def __add__(self, other):
        if (self.xdim, self.ydim) != (other.xdim, other.ydim):
            raise ValueError("Addition of matrices with dimensions {(self.xdim,self.ydim)} and {(other.xdim,other.ydim)} not defined.")

        return SymbolicMatrix(tuple(self.matrix[i][j]+other.matrix[i][j] for i in range(self.xdim)) for j in range(self.ydim))

    def __sub__(self, other):
        if (self.xdim, self.ydim) != (other.xdim, other.ydim):
            raise ValueError("Subtraction of matrices with dimensions {(self.xdim,self.ydim)} and {(other.xdim,other.ydim)} not defined.")

        return SymbolicMatrix(tuple(self.matrix[i][j]-other.matrix[i][j] for i in range(self.xdim)) for j in range(self.ydim))

    def __mul__(self, other):
        return SymbolicMatrix(tuple(sum(self.matrix[i][k]*other.matrix[k][j] for k in range(self.ydim)) for j in range(other.ydim)) for i in range(self.xdim))

    def __pow__(self, it):
        if it < 0 or not isinstance(it,int):
            raise ValueError("Power must be integer >= 0")
        elif self.xdim != self.ydim:
            raise ValueError("Matrix exponent only defined for square matrices.")
        elif it == 0:
            return self.identity()
        elif it == 1:
            return SymbolicMatrix(self.matrix)
        elif it == 2:
            return self*self
        # use divide and conquer algorithm for larger powers
        else:
            return self**(-(-it // 2)) * self**(it // 2)

    def __str__(self):
        max_w = [max(len(str(sublist[i])) for sublist in self.matrix) for i in range(len(self.matrix[0]))]
        return "[" + ",\n ".join("["+', '.join(f"{str(item):{max_w[i]}}" for i,item in enumerate(sublist))+"]" for sublist in self.matrix) + "]"

    #  def as_latex(self):
        # TODO: finish this function, and write methods in .numerics as well
        #  str1 = r"\begin{pmatrix}"
        #  str2 = "\n".join("   " + "&".join(s.as_latex() for s in sublist) + r"\\" for sublist in self.matrix)
