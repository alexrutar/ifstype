import typing
import operator
from sympy import Rational

# TODO: implement arbitrary unions / intersections of intervals, have a multiple interval object, etc.
class IntervalSet:
    def __init__(self, *intervals):
        # process intervals
        self.interval = interval

# in order to suppose NetInterval objects, Intervals have an additional _param attribute that is never observed from outside
# TODO: this is garbage, fix it
class Interval(tuple):
    "A fast interval class"
    a = property(operator.itemgetter(0))
    b = property(operator.itemgetter(1))
    def __new__(cls,a: typing.Optional[Rational] = None, b: typing.Optional[Rational] = None, _param: typing.Optional[Rational] = Rational(0)):
        if a is None or b is None or a > b:
            return super().__new__(cls,(None,None,_param))
        else:
            return super().__new__(cls,(a,b,_param))

    def copy(self):
        return Interval(self.a,self.b)

    def is_empty(self):
        return self.a is None and self.b is None

    def delta(self):
        if self.is_empty():
            return 0
        else:
            return self.b-self.a

    def is_point(self):
        return not self.is_empty() and self.a == self.b

    # contains
    def __contains__(self, item: Rational):
        return self.a <= item and self.b >= item
    def subset(self, other):
        return (not other.is_empty()) and self.a >= other.a and other.b >= self.b
    def supset(self, other):
        return other.is_empty() or (self.a <= other.a and other.b <= self.b)
    def proper_subset(self,other):
        return self.subset(other) and self != other
    def proper_supset(self,other):
        return self.supset(other) and self != other

    # representation
    def __str__(self):
        return f"Iv[{self.a},{self.b}]"

    def __repr__(self):
        return f"Iv[{self.a},{self.b}]"

    def __mul__(self, cnst: Rational):
        if self.is_empty():
            return Interval()
        elif cnst < 0:
            return Interval(self.b * cnst, self.a * cnst)
        else:
            return Interval(self.a * cnst, self.b * cnst)

    def __truediv__(self,cnst: Rational):
        if self.is_empty():
            return Interval()
        elif cnst < 0:
            return Interval(self.b / cnst, self.a / cnst)
        else:
            return Interval(self.a / cnst, self.b / cnst)

    def __add__(self,cnst: Rational):
        if self.is_empty():
            return Interval()
        else:
            return Interval(self.a + cnst, self.b + cnst)

    def __sub__(self,cnst: Rational):
        if self.is_empty():
            return Interval()
        else:
            return Interval(self.a - cnst, self.b - cnst)

    def __and__(self,other: 'Interval'):
        "intersection"
        if self.is_empty() or other.is_empty():
            return Interval()
        else:
            return Interval(max(self.a,other.a),min(self.b,other.b))

    # caution! cannot support union
    def open_intersect(self,other):
        intersection = self & other
        return not(intersection.is_empty() or intersection.is_point())


class NetInterval(Interval):
    "A special interval type representing a net interval of generation alpha"
    # TODO: fix persistent mutability of alpha
    alpha = property(operator.itemgetter(2))
    def __new__(cls,alpha,a,b):
        self = super().__new__(cls,a,b,alpha)
        return self

    # representation
    def __str__(self):
        return "NetIv({})[{},{}]".format(self.alpha,self.a,self.b)
    def __repr__(self):
        return "NetIv({})[{},{}]".format(self.alpha,self.a,self.b)


