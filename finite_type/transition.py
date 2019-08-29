import graph_tool as gt
import numpy as np
import typing
from functools import reduce
import itertools
from operator import mul

from .numerics.rational import Constants as C, Rational
from .ifs import CtrFunc

class Neighbour(CtrFunc):
    @classmethod
    def from_f(cls,f,iv):
        func = CtrFunc(iv.delta,iv.a).inverse().compose(f)
        return cls(func.r,func.d)

    def to_f(cls, net_iv):
        return net_iv.normalization_func().compose(self)

    @property
    def L(self):
        return self.r

    @property
    def a(self):
        return -self.d

class NeighbourSet(tuple):
    """
    A neigbour set is just a sorted tuple of neighbours (contraction functions).
    """
    def __new__(cls,nb_itbl):
        self = super().__new__(cls,sorted(set(nb_itbl)))
        self.lmax = max(abs(nb.L) for nb in self)
        return self

    def __str__(self):
        return ", ".join(f"({nb.d},{nb.L})" for nb in self)

    @classmethod
    def base(cls):
        return cls((Neighbour(d=0,r=1),))

    def maximal_nbs(self):
        """Compute the neigbours in self.nb_set of maximal size"""
        return (nb for nb in self if abs(nb.L) == self.lmax)

    def nonmaximal_nbs(self):
        """Compute the neigbours in self.nb_set not of maximal size"""
        return (nb for nb in self if abs(nb.L) != self.lmax)

class TransitionMatrix:
    """TransitionMatrix class to represent the transition matrix associated to an edge in the transition graph.
    """
    def __init__(self, double_list,exact=False):
        self.exact = exact
        if exact:
            self.matrix = np.array(double_list,dtype="object")
            self.matrix.flags.writeable=False
        else:
            self.matrix = np.array(double_list,dtype=float)

    def pos_row(self):
        return all(any(x>0 for x in row) for row in self.matrix)

    def spectral_radius(self):
        return np.abs(np.linalg.eigvals(self.matrix.astype(float))).max()

    def __hash__(self):
        return hash(self.matrix.tobytes())

    def __repr__(self):
        return repr(self.matrix)

    def __mul__(self, other):
        return TransitionMatrix(self.matrix.dot(other.matrix),self.exact and other.exact)

    def __repr__(self):
        return repr(self.matrix)

    def __str__(self):
        max_w = [max(len(str(sublist[i])) for sublist in self.matrix) for i in range(len(self.matrix[0]))]
        return "[" + ",\n ".join("["+','.join(f"{str(item):{max_w[i]}}" for i,item in enumerate(sublist))+"]" for sublist in self.matrix) + "]"

    #  def as_latex(self):
        # TODO: finish this function, and write methods in .numerics as well
        #  str1 = r"\begin{pmatrix}"
        #  str2 = "\n".join("   " + "&".join(s.as_latex() for s in sublist) + r"\\" for sublist in self.matrix)

class EdgeInfo(typing.NamedTuple):
    """Named tuple to represent the edge information intrinsic to an edge (other than the source and target)"""
    t_index: Rational # if speed is ever an issue, this can be converted to an index int
    length: Rational
    transition: TransitionMatrix

class LocalDim:
    """Class to represent local dimensions, for easy printing / visualization of what the values are."""
    def __init__(self, spr, ln):
        self.spr = spr
        self.length = ln

    def __float__(self):
        return float(np.log(float(self.spr))/np.log(float(self.length)))

    def __str__(self):
        return f"log({self.spr})/log({self.length})"

class TransitionGraph:
    def __init__(self):
        self.is_wft = False
        self.g = gt.Graph()

        # graph properties
        self._edge_info = self.g.new_edge_property("object")
        self.nb_set_lookup = self.g.new_vertex_property("object") # correspondence of vertex objects to neighbour sets

        # additional associations
        self.vtx_lookup = {} # correspondence of neighbour set objects to vertices
        self.edge_lookup = {} # correspondence of edge indices to edge descriptors

    def all_nb_sets(self):
        "Iterable for all neighbour sets, which also allows inclusion checking."
        return self.vtx_lookup.keys()

    # -------------------------------------------------
    # graph information methods
    # -------------------------------------------------
    def get_identifier(self, nb_set):
        "Get a unique identifier corresponding to the neighbour set."
        return str(self.vtx_lookup[nb_set])

    def get_nb_set(self, vtx):
        "Get the neighbour set corresponding to a given vertex."
        return self.nb_set_lookup[vtx]

    def has_pos_row(self):
        "Check if the transition graph has the positive row property; in other words, that ever transition matrix has a positive entry in every row."
        return all(e_inf.transition.pos_row() for e_inf in self._edge_info)

    def edge_info(self,e,by_label=False):
        """Return the EdgeInfo object associated with the edge."""
        if by_label:
            return self._edge_info[self.edge_lookup[e]]
        else:
            return self._edge_info[e]

    def edge_info_string(self, e, by_label=False):
        "Return an information string for attributes of the given edge."
        return f"{e.source()} -> {e.target()}\nindex={self.edge_info(e,by_label=by_label).t_index}\nlength={self.edge_info(e).length}\nlabel={self.g.edge_index[e]}\n" + f"{self.edge_info(e).transition}\n"


    # -------------------------------------------------
    # compute local dimension attributes
    # -------------------------------------------------
    def is_valid_path(self, path, by_label=False):
        """Check if a sequences if edges is in fact a valid path."""
        return all(e1.target() == e2.source() for e1,e2 in zip(path,path[1:]))

    def local_dim(self, loop, by_label=True):
        """Compute the local dimension associated with a loop."""
        if by_label:
            loop = [self.edge_lookup[label] for label in loop]
            assert self.is_valid_path(loop), f"{loop} is not a valid path"
            assert loop[-1].target() == loop[0].source(), f"{loop} start vertex and end vertex are distinct"

        mat = reduce(mul,(self._edge_info[e].transition for e in loop))
        L = reduce(mul, (self._edge_info[e].length for e in loop))

        return LocalDim(mat.spectral_radius(),L)

    def vertex_local_dims(self,vtx,search_depth=1):
        """Returns a sorted list of all possible local dimensions for a cycle (no repeated vertices) containing vtx."""
        cycles = gt.topology.all_paths(self.g,vtx,vtx,edges=True)
        cycle_prods = itertools.product(cycles, repeat=search_depth)
        return sorted(float(self.local_dim(list(itertools.chain.from_iterable(loop_tups)),by_label=False)) for loop_tups in cycle_prods)
    

    # -------------------------------------------------
    # compute topological aspects of the graph
    # -------------------------------------------------
    def essential_class(self):
        """Returns a graph_view consisting of the essential class of the graph."""
        components, _, att = gt.topology.label_components(self.g,attractors=True) # get the components and the attractor
        attractor = np.array([(att == True)[comp] == 1 for comp in components.a])
        return gt.GraphView(self.g,vfilt=attractor)

    def components(self):
        """Create a vertex property map `vprop` such that for any v,
        - vprop[v] == -1 if v is an isolated component with no self loop
        - vprop[v] == 0 if v is in the attractor (essential class if weak finite type)
        - vprop[b] == n for 1 <= n <= N where N is the number of non-trivial components of the graph.
        """
        components, hist, att = gt.topology.label_components(self.g,attractors=True) # get the components and the attractor
        
        # boolean array for attractor of the graph
        attractor = np.array([(att == True)[comp] == 1 for comp in components.a])

        # boolean array for isolated vertices of the graph
        isolated= np.array([(hist == 1)[components[vtx]] and not vtx in vtx.out_neighbours() for vtx in self.g.vertices()],dtype=bool)

        # set isolated components to be type -1
        np.place(components.a, isolated, -1)
        
        # set attractor to be type 0
        np.place(components.a, attractor, 0)

        # set all other loop classes to be type 1 through n
        uniques = np.unique(components.a)
        lookup_map = {elem:idx for idx,elem in enumerate(uniques[uniques > 0],1)}
        components.a = np.array([lookup_map.get(x,x) for x in components.a])

        return components


    # -------------------------------------------------
    # internal construction methods
    # -------------------------------------------------
    def add_nb_set(self, nb_set, transition=None, t_index=None, prec=None,length=None):
        if nb_set not in self.vtx_lookup.keys():
            v = self.g.add_vertex()
            self.nb_set_lookup[v] = nb_set
            self.vtx_lookup[nb_set] = v
        if prec in self.all_nb_sets():
            e = self.g.add_edge(self.vtx_lookup[prec],self.vtx_lookup[nb_set])
            self.edge_lookup[self.g.edge_index[e]] = e # register edge by index
            self._edge_info[e] = EdgeInfo(
                    t_index=t_index,
                    transition=transition,
                    length=length)

    def remove_terminal_vertices(self):
        """Remove vertices with out degree 0, repair the vtx_lookup, and return the corresponding list of neighbour sets removed.
        Warning: this changes edge numbering."""
        # TODO: use graph views / filters instead of mutating the graph
        zero_degs = gt.util.find_vertex(self.g,"out",0)
        zdg = set()
        while(len(zero_degs) > 0):
            zdg.update(self.nb_set_lookup[v] for v in zero_degs)
            self.g.remove_vertex(zero_degs)
            zero_degs = gt.util.find_vertex(self.g,"out",0)

        # rebuild vertex lookup
        self.vtx_lookup = {self.nb_set_lookup[v]:v for v in self.g.vertices()}
        return zdg
