
class Interval:
    def __init__(self,a,b):
        if a>b:
            self.left = None
            self.right = None
            self.is_empty = True
            self.is_point = False
        else:
            self.left = a
            self.right = b
            self.is_empty = False
            if a == b:
                self.is_point = True
            else:
                self.is_point = False

    @classmethod
    def empty(cls):
        return cls(-1,0)

    def as_latex(self):
        return "{}/{}".format(float(self.left),float(self.right))

    def __str__(self):
        return str((self.left,self.right))

    def __repr__(self):
        return str((self.left,self.right))

    def __eq__(self, other):
        return self.left == other.left and self.right == other.right

    def __hash__(self):
        return hash((self.left,self.right))

    def __mul__(self, cnst):
        if self.is_empty:
            return Interval.empty()
        elif cnst < 0:
            return Interval(self.right * cnst, self.left * cnst)
        else:
            return Interval(self.left * cnst, self.right * cnst)

    def __truediv__(self,cnst):
        if self.is_empty:
            return Interval.empty()
        elif cnst < 0:
            return Interval(self.right / cnst, self.left / cnst)
        else:
            return Interval(self.left / cnst, self.right / cnst)
        
    def __add__(self,cnst):
        if self.is_empty():
            return Interval.empty()
        else:
            return Interval(self.end[0] + cnst, self.end[1] + cnst)

    def __sub__(self,cnst):
        if self.is_empty():
            return Interval.empty()
        else:
            return Interval(self.end[0] - cnst, self.end[1] - cnst)

    def __and__(self,other):
        "Set intersection"
        if self.is_empty or other.is_empty:
            return Interval.empty()
        elif self.left <= other.left:
            return Interval(other.left,self.right)
        else:
            return Interval(self.left,other.right)

    def __or__(self,other):
        if self.is_empty or other.is_empty:
            return Interval.empty
        else:
            return Interval(min(self.left,other.left),max(self.right,other.right))
    
    def len(self):
        if self.is_empty or self.is_point:
            return 0
        else:
            return self.right-self.left


def open_overlap(i1,i2):
    intersection = i1.intersect(i2)
    return not(intersection.is_empty() or intersection.is_point())

def neighbours(interval,level_set):
    return {iv for iv in level_set if open_overlap(interval,iv)}

def normal_neighbours(interval,level_set):
    return {(iv-interval.end[0])*(1/interval.meas()) for iv in level_set if open_overlap(interval,iv)}

#  def nbhds(interval_set):
    #  return eq_part(interval_set, open_overlap)

#  def eq_part(iterable, relation):
    #  """Partitions a set of objects into equivalence classes

    #  Args:
        #  iterable: collection of objects to be partitioned
        #  relation: equivalence relation. I.e. relation(o1,o2) evaluates to True
            #  if and only if o1 and o2 are equivalent

    #  Returns: classes, partitions
        #  classes: A sequence of sets. Each one is an equivalence class
        #  partitions: A dictionary mapping objects to equivalence classes
    #  """
    #  classes = []
    #  partitions = {}
    #  for o in iterable:
        #  # find the class it is in
        #  found = False
        #  for c in classes:
            #  if relation(next(iter(c)), o):  # is it equivalent to this class?
                #  c.add(o)
                #  partitions[o] = c
                #  found = True
                #  break
        #  if not found:  # it is in a new class
            #  classes.append(set([o]))
            #  partitions[o] = classes[-1]
    #  return classes, partitions

