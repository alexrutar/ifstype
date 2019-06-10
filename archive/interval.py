import itertools


class Word:
    # must have length at least 1
    def __init__(self, ab, r, rm):
        self.r = r
        self.rm = rm
        self.ab = tuple(ab)
        self.len = len(ab)

    def __hash__(self):
        return hash(self.ab)

    def __str__(self):
        return str(self.ab)

    def __repr__(self):
        return str((self.ab,self.r,self.rm))

    @classmethod
    def empty(cls):
        return cls("",1,None)

    def extend(self, a, rnew):
        new_r = self.r * rnew
        new_ab = self.ab + (a,)
        if self.len == 0:
            new_rm = 1
        else:
            new_rm = self.r
        return Word(new_ab, new_r,new_rm)


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
        return cls(-1,0)

    # contains
    def __contains__(self, item):
        return self.a <= item and self.b >= item
    # representation
    def __str__(self):
        return str((self.a,self.b))

    def __repr__(self):
        return str((self.a,self.b))

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
        if self.is_empty or other.is_empty:
            return Interval.empty
        else:
            return Interval(min(self.a,other.a),max(self.b,other.b))

def draw_level(idx,ivs):
   return "        \\foreach \\x/\\y in {" + ",".join([w.as_latex() for w in ivs]) + "} {\n" +\
          "            \\draw (\\x,{}) -- (\\y,{});\n".format(idx,idx) +\
          "            \\draw (\\x,{}-0.2) -- (\\x,{}+0.2);\n".format(idx,idx) +\
          "            \\draw (\\y,{}-0.2) -- (\\y,{}+0.2);\n".format(idx,idx) +\
          "        }\n"


from collections import defaultdict
class NetIntervals:
    def __init__(self,words,ifs):
        # takes list((iv,r),...)
        # keep track of intervals, and to each interval, associate the words which generate it (can be many)
        self.words = set(words)
        unique = sorted(((ifs.word_interval(w), w) for w in self.words),key=lambda x:x[0])
        endpoints = itertools.chain.from_iterable([(iv.a,w),(iv.b,w)] for iv,w in unique)

        # creates a list of intervals and a dictionary showing which words generate a given interval
        self.iv_dict = defaultdict(set)
        for k,v in unique:
            self.iv_dict[k].add(v)
        self.iv = sorted(self.iv_dict.keys())

        # creates a list of endpoints and a dictionary showing which words generate a given endpoint
        self.end_dict = defaultdict(set)
        for k,v in endpoints:
            self.end_dict[k].add(v)
        self.end = sorted(self.end_dict.keys())

        # the net intervals from the words
        self.net = [Interval(a,b) for a,b in zip(self.end,self.end[1:])]

        # all the possible neighbourhoods and an assignment of type numbers to the neighbourhoods
        self.nbhd = {iv:self.nbhd(iv) for iv in self.net}
        self.nbhd_types={nb:idx for idx,nb in enumerate(set(self.nbhd.values()))}

    def min_delta(self):
        return min(iv.delta for iv in self.net)

    def max_delta(self):
        return max(iv.delta for iv in self.net)

    def nbhd(self,interval,normalize_factor=None):
        "Compute the neighbour set of a given interval"
        if interval.is_empty or interval.is_point:
            return []
        if normalize_factor is None:
            nf = interval.delta
        else:
            nf = normalize_factor
        return tuple(sorted((iv-interval.a)/nf for iv in self.iv if open_overlap(iv,interval)))

    # old
    #  def nbhd_types(self,normalize_factor=None):
        #  "Generate all neighbourhood types"
        #  # also generates a dictionary mapping net intervals to neighbourhood types
        #  if normalize_factor is None:
            #  f = lambda x: x.delta
        #  else:
            #  f = lambda x: normalize_factor
        #  self.nbhd = {iv:self.nbhd(iv,normalize_factor=f(iv)) for iv in self.net}
        #  return set(self.nbhd.values())

    def print_as_latex(self):
        levels = [[]]
        for iv in self.iv:
            done = False
            for level in levels:
                if level == [] or level[-1].b <= iv.a:
                    done = True
                    level.append(iv)
                    break
            if not done:
                levels.append([iv])

        print("\\begin{center}\n" + 
              "    \\begin{tikzpicture}[xscale=14,yscale=0.2]\n")
        for idx, level in enumerate(levels):
            print(draw_level(idx,level))
        print("    \\end{tikzpicture}\n"
              "\\end{center}\n")


def open_overlap(i1,i2):
    intersection = i1 & i2
    return not(intersection.is_empty or intersection.is_point)

