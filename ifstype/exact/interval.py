import typing
import attr
import numbers

from .rational import Constants as C


# must check a <= b
@attr.s(frozen=True,slots=True)
class Interval:
    "A fast closed interval class"
    a: numbers.Real = attr.ib(default=C.n_0)
    b: numbers.Real = attr.ib(default=C.n_1)

    # -----------------------------------------------
    # construction methods
    # -----------------------------------------------
    @classmethod
    def point(cls,a=None):
        return cls(a,a)

    def copy(self):
        """Create a new instance of the same interval type but with new endpoints"""
        return Interval(a=a,b=b)
    
    # -----------------------------------------------
    # standard interval properties
    # -----------------------------------------------
    @property
    def is_empty(self):
        return self.a > self.b

    @property
    def is_point(self):
        return self.a == self.b

    @property
    def delta(self):
        if self.is_empty:
            return 0
        else:
            return self.b-self.a

    # -----------------------------------------------
    # inclusion methods
    # -----------------------------------------------
    def __contains__(self, item):
        return (not self.is_empty) and a<= item and item <= b
    def subset(self, other):
        return self.is_empty or ((not other.is_empty) and self.a >= other.a and other.b >= self.b)
    def supset(self, other):
        return other.subset(self)
    def proper_subset(self,other):
        return self.subset(other) and self != other
    def proper_supset(self,other):
        return other.subset(self) and self != other


    # -----------------------------------------------
    # misc methods
    # -----------------------------------------------
    def __str__(self):
        if self.is_empty:
            return "EmptyInterval"
        elif self.is_point:
            return f"PointInterval({self.a})"
        else:
            return f"Interval(a={self.a},b={self.b})"

    def __repr__(self):
        return f"Interval(a={self.a},b={self.b})"

