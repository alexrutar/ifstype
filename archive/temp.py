class Interval:

    @classmethod
    def empty(cls):
        return cls(0,-1,None)

    def as_latex(self):
        return "{}/{}".format(float(self.a),float(self.b))

    # representation
    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        if self.is_empty:
            return "()"
        else:
            return "{}{},{}{}".format(self.t[0],self.a,self.b,self.t[1])

    def __hash__(self):
        return hash((self.a,self.b,self.t))

def iv_empty(n):
    def iv_internal(f):
        def wrapper(*args):
            if any(args[i].is_empty for i in range(n)):
                return Interval.empty()
            else:
                return f(*args)
        return wrapper
    return iv_internal

class IntervalCC:
    def __init__(self,a,b):
        if a>b:
            self.a = None
            self.b = None
            self.is_empty = True
        else:
            self.a = a
            self.b = b
            self.is_empty = False

    # contains
    def __contains__(self, item):
        return self.a <= item and self.b >= item

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

    @iv_empty(1)
    def __mul__(self, cnst):
        if cnst < 0:
            return IntervalCC(self.b * cnst, self.a * cnst)
        else:
            return IntervalCC(self.a * cnst, self.b * cnst)

    @iv_empty(1)
    def __truediv__(self,cnst):
        if cnst < 0:
            return IntervalCC(self.b / cnst, self.a / cnst)
        else:
            return IntervalCC(self.a / cnst, self.b / cnst)

    @iv_empty(1)
    def __add__(self,cnst):
        return IntervalCC(self.a + cnst, self.b + cnst)

    @iv_empty(1)
    def __sub__(self,cnst):
        else:
            return IntervalCC(self.a - cnst, self.b - cnst)

    @iv_empty(2)
    def __and__(self,other):
        "intersection"
        if self.a <= other.a:
            return IntervalCC(other.a,self.b)
        else:
            return IntervalCC(self.a,other.b)
    @iv_empty(2)
    def __or__(self,other):
        "union"
        return IntervalCC(min(self.a,other.a),max(self.b,other.b))


