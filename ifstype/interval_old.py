import typing
from sympy import Rational

# TODO: implement arbitrary unions / intersections of intervals, have a multiple interval object, etc.
class IntervalSet:
    def __init__(self, *intervals):
        # process intervals
        self.interval = interval

class Interval:
    def __init__(self,a: typing.Optional[Rational], b: typing.Optional[Rational]):
        if a>b:
            self.a = None
            self.b = None
            self.is_empty = True
            self.is_point = False
            self.delta = 0
        else:
            self.a = a
            self.b = b
            self.is_empty = False
            if a == b:
                self.is_point = True
                self.delta = 0
            else:
                self.is_point = False
                self.delta = b-a

    def as_latex(self):
        return "{}/{}".format(float(self.a),float(self.b))

    def ep(self):
        return (self.a,self.b)
    @classmethod
    def empty(cls):
        return cls(0,-1)

    # contains
    def __contains__(self, item: Rational):
        return self.a <= item and self.b >= item
    def subset(self, other):
        return self.a >= other.a and other.b >= self.b
    def supset(self, other):
        return self.a <= other.a and other.b <= self.b

    # representation
    def __str__(self):
        return "Iv[{},{}]".format(self.a,self.b)

    def __repr__(self):
        return "Iv[{},{}]".format(self.a,self.b)

    # math properties
    def __eq__(self, other):
        return self.a == other.a and self.b == other.b
    def __lt__(self, other):
        return self.a < other.a or (self.a == other.a and self.b < other.b)
    def __gt__(self, other):
        return self.a > other.a or (self.a == other.a and self.b > other.b)
    def __le__(self, other):
        return self.a <= other.a or (self.a == other.a and self.b <= other.b)
    def __ge__(self, other):
        return self.a >= other.a or (self.a == other.a and self.b >= other.b)

    def __hash__(self):
        return hash((self.a,self.b))

    def __mul__(self, cnst):
        if self.is_empty:
            return Interval.empty()
        elif cnst < 0:
            return Interval(self.b * cnst, self.a * cnst)
        else:
            return Interval(self.a * cnst, self.b * cnst)
    def __truediv__(self,cnst):
        if self.is_empty:
            return Interval.empty()
        elif cnst < 0:
            return Interval(self.b / cnst, self.a / cnst)
        else:
            return Interval(self.a / cnst, self.b / cnst)
    def __add__(self,cnst):
        if self.is_empty:
            return Interval.empty()
        else:
            return Interval(self.a + cnst, self.b + cnst)
    def __sub__(self,cnst):
        if self.is_empty:
            return Interval.empty()
        else:
            return Interval(self.a - cnst, self.b - cnst)
    def __and__(self,other):
        "intersection"
        if self.is_empty or other.is_empty:
            return Interval.empty()
        elif self.a <= other.a:
            return Interval(other.a,self.b)
        else:
            return Interval(self.a,other.b)
    def __or__(self,other):
        "union"
        if self.is_empty and other.is_empty:
            return Interval.empty()
        elif self.is_empty:
            return other
        elif other.is_empty:
            return self
        else:
            return Interval(min(self.a,other.a),max(self.b,other.b))

    def open_intersect(self,other):
        intersection = self.__and__(other)
        return not(intersection.is_empty or intersection.is_point)


class NetInterval(Interval):
    "A special interval type representing a net interval of generation alpha"
    def __init__(self,alpha,a,b):
        super().__init__(a,b)
        self.alpha = alpha
    # representation
    def __str__(self):
        return "NetIv({})[{},{}]".format(self.alpha,self.a,self.b)

    def __repr__(self):
        return "NetIv({})[{},{}]".format(self.alpha,self.a,self.b)


