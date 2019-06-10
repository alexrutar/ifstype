from collections import defaultdict
import itertools
import random

class Interval:
    def __init__(self,a,b):
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

    @classmethod
    def empty(cls):
        return cls(0,-1)

    # contains
    def __contains__(self, item):
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


class IntervalNet:
    """An IntervalNet takes a generating set of words and an optional range on which the net intervals are valid
    Assumes that elements of 'words' all generate Interval(0,1)
    """
    def __init__(self,alpha,words,interval=Interval(0,1)):
        # warning: open intersection is not sufficient, since contributions to endpoints is also important!
        self.words = set(w for w in words if not (w.interval() & interval).is_empty) # only keep words which intersect with interval
        self.alpha = alpha
        self.interval = interval

        # creates a list of intervals and a dictionary showing which words generate a given interval
        self.iv_dict = defaultdict(set)
        for w in self.words:
            self.iv_dict[w.interval()].add(w)
        self.iv = sorted(self.iv_dict.keys())

        #  print(sorted(set(itertools.chain.from_iterable([iv.a,iv.b] for iv in self.iv))))
        endpoints = sorted(ep for ep in set(itertools.chain.from_iterable([iv.a,iv.b] for iv in self.iv)) if interval.a <= ep and ep <= interval.b)

        # the net intervals from the words
        self.net = sorted(set(NetInterval(alpha,a,b) for a,b in zip(endpoints,endpoints[1:])))

        # all the possible neighbour sets
        self.nb_set_dict= {net_iv:tuple(sorted(self.norm(iv,net_iv) for iv in self.iv if iv.open_intersect(net_iv))) for net_iv in self.net}
        self.nb_set_types = set(self.nb_set_dict.values())

    def __iter__(self):
        return iter(self.net)

    def adjacent(self, net_iv):
        "Compute the adjacent intervals (left,right) if they exist (returns empty if they don't)"
        idx = self.net.index(net_iv)
        try:
            left = self.net[idx-1]
        except IndexError:
            left = Interval.empty()
        try:
            right = self.net[idx+1]
        except IndexError:
            right = Interval.empty()
        return (left,right)

    def random_net_iv(self,n=1):
        "Get a random net interval"
        return random.sample(self.net, n)

    def norm(self,iv,net_iv):
        "normalize the interval intersecting with net_iv"
        return (iv-net_iv.a)/net_iv.delta

    def norminv(self,nb_iv,net_iv):
        "inverse normalize an element of the neighbour set of net_iv to get the corresponding interval"
        return nb_iv*net_iv.delta+net_iv.a

    def iv_words(self,iv):
        "Get the the words which generate a given interval"
        return self.iv_dict[iv]

    def nb_set(self,net_iv):
        "Return the local neighbour set of a net_iv"
        return self.nb_set_dict[net_iv]

    def argmin_delta(self):
        "Returns the net interval with minimum size"
        return min(self.net, key=lambda x:x.delta)

    def word_contributions(self,net_iv):
        "Returns an iterable of (interval,word) pairs which contribute to net_iv (for recursive calculation purposes)"
        return ((iv,word) for iv in self.nb_set(net_iv) for word in self.iv_words(self.norminv(iv,net_iv)))

