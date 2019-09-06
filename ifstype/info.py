import matplotlib.pyplot as plt
import matplotlib
import graph_tool as gt
import numpy as np
import itertools
from collections import defaultdict

from .exact.algebraic import AlgebraicNumber

class GraphWriter:
    """Class to create human readable information strings of the TransitionGraph and the associated iterated function system."""
    def __init__(self, tr_graph):
        self.tr_graph = tr_graph

    def variable_info(self):
        all_nums = itertools.chain((f.r for f in self.tr_graph.ifs.funcs),(f.d for f in self.tr_graph.ifs.funcs))
        # pull a number field, if it exists
        num_field = None
        for n in all_nums:
            if isinstance(n,AlgebraicNumber):
                num_field = n.num_field

        if num_field is not None:
            return f"'{num_field.symbol}' is an algebraic number satisfying '{num_field.minpoly.with_symbol(num_field.symbol)} = 0' with approximate value '{num_field.expr_float}'."
        else:
            return ""


    def ifs_info(self):
        def as_str(num):
            if isinstance(num,AlgebraicNumber):
                return num.as_str(term=True,space=False)
            else:
                return str(num)
                
        return "- " + "\n- ".join(f"S{i}(x)=x*{as_str(f.r)} + {as_str(f.d)}" for i,f in enumerate(self.tr_graph.ifs.funcs,1))

    def nb_set_info(self):
        return "\n".join(f"{self.tr_graph.get_identifier(nb_set)} : {nb_set}" for nb_set in self.tr_graph.all_nb_sets())

    def edge_info(self, e, by_label=True):
        return f"{e.source()} -> {e.target()}\nedge index={self.tr_graph.edge_info(e,by_label=by_label).t_index}\nlength={self.tr_graph.edge_info(e).length}\nlabel={self.tr_graph.g.edge_index[e]}\n" + f"{self.tr_graph.edge_info(e).transition}\n"

    def all_edge_info(self):
        return "\n".join(self.edge_info(e, by_label=False) for e in self.tr_graph.g.edges())

    def info(self):
        return "".join(header + ":\n" + content + "\n\n\n"for header, content in [("Variables",self.variable_info()), ("Contraction Functions", self.ifs_info()), ("Neighbour Sets",self.nb_set_info()), ("Edge Information", self.all_edge_info())])

    def info_to_file(self,filename):
        # write the neighbour sets
        with open(filename,'w') as outfile:
            outfile.write(self.info())



class GraphArtist:
    """Class to create .pdf visualizations of the TransitionGraph object."""
    def __init__(self,tr_graph):
        self._tr_graph = tr_graph
        self.colours = {
                "base_light":'#E9EDED',
                "base_dark":'gray',
                "strong_highlight" : '#C40E0E',
                "blue" : 'blue',
                "green" : 'green',
                "black" : 'black'}

        self.base_light_color = '#E9EDED'
        self.highlight_color = '#C40E0E'
        self.blue_color = 'blue'
        self.cmap = matplotlib.cm.get_cmap('viridis')
        self.scale=1

    # methods to draw the interval class
    def _resolve_level(self, iv_itbl):
        """given an iterable of intervals, _resolve_levels creates a sorted list of the unique items, and then places them in sublist.
        Within each sublist, no two intervals overlap (except possibly at an endpoint).
        The number of sublists is the maximum number of overlaps of interiors of intervals.
        """
        level = [[]]
        for iv in sorted(set(iv_itbl)):
            done = False
            for layer in level:
                if layer == [] or layer[-1].b <= iv.a:
                    done = True
                    layer.append(iv)
                    break
            if not done:
                level.append([iv])
        return level

    def _set_scale(scale, net_iv_set):
        min_width = min(net_iv.delta for net_iv in net_iv_set)
        # bounds
        max_scale = 18
        min_scale = 1

        if isinstance(scale,(float,int)):
            self.scale = scale
        elif scale == "relative":
            self.scale = .01/float(min_width)
        elif scale == "absolute":
            self.scale = 1
        elif scale == "wide":
            self.scale = .015/float(min_width)
        elif scale == "maximum":
            self.scale = max_scale
        else:
            raise ValueError("Invalid scale argument.")

        # correction within bounds
        if self.scale > max_scale:
            self.scale = max_scale
        elif self.scale < min_scale:
            self.scale = min_scale

    def net_intervals(self, net_iv_set,filename,scale="relative"):
        """Draw the generations image of the interval Generations object.
        The optional argument 'scale' denotes how large to make the image:
        - relative: increase size to exactly fit the smallest net interval
        - wide: fit the smallest net interval with some extra space
        - absolute: always the same (not too large) size
        """
        # compute endpoint locations, and save net intervals
        all_eps = sorted(set(itertools.chain((net_iv.a for net_iv in net_iv_set),(net_iv.b for net_iv in net_iv_set))))
        ep_locs = [float(s) for s in all_eps]
        ep_labels = [str(s) for s in all_eps]

        levels = defaultdict(set)
        for net_iv in net_iv_set:
            levels[net_iv.alpha].add(net_iv)

        # TODO: draw negative intervals with a different colour
        # resolve interval levels
        ivl_draws = {alpha:self._resolve_level(itertools.chain.from_iterable((f.interval() for f in net_iv.containing_funcs()) for net_iv in net_ivs_layer)) for alpha,net_ivs_layer in levels.items()}
        heights = {alpha:len(level) for alpha, level in ivl_draws.items()}

        # compute dimensions
        step_height = 0.35

        net_place = {}
        total_height = 0
        for alpha in sorted(heights.keys()):
            net_place[alpha] = total_height
            total_height += (heights[alpha]+1.8)* step_height

        xgap = 0.05/self.scale
        ydims = (-step_height,total_height+step_height)
        xdims = (-xgap,1+xgap)

        # figure creation
        fig = plt.figure(figsize=(self.scale*36.0,2*total_height))
        ax = fig.add_subplot()
        ax.set_axisbelow(True)

        # set vertical dashed line indicators and x-labels
        plt.vlines(ep_locs,*ydims,ls=':',color=self.colours["base_light"])
        plt.xticks(ep_locs,ep_labels,rotation="vertical")

        # set alpha labels and dashed line separators
        for alpha, ht in net_place.items():
            plt.text(-0.04/self.scale,ht,f"α = {alpha}",ha='left',va='baseline')
            plt.hlines(ht-0.2,*xdims,ls=':',color=self.colours["base_dark"])

        for alpha, level in ivl_draws.items():
            # plot the net intervals
            for net_iv in levels[alpha]:
                self._draw_net_iv(net_place[alpha], net_iv)
            # plot the intervals containing them
            for idx,layer in enumerate(level):
                for ivl in layer:
                    self._draw_interval(net_place[alpha] + (1+idx)*step_height, ivl, color=self.colours["black"])


        # no y-axis
        ax.get_yaxis().set_visible(False)

        # set dimensions
        plt.xlim(*xdims)
        plt.ylim(*ydims)

        # save figure with tight boundary
        plt.savefig(filename, bbox_inches = 'tight')

    def _draw_interval(self, ht, iv, color='black',mid_label=None):
        """Plot an interval with y value ht."""
        lw = 2
        iv_height = 0.1
        plt.hlines(ht, iv.a, iv.b, color, lw=lw)
        plt.vlines([iv.a,iv.b], ht+iv_height, ht-iv_height, color, lw=lw)
        if mid_label is not None:
            plt.text((iv.a+iv.b)/2,ht,mid_label,ha='center',va='center',bbox={'facecolor':'white','edgecolor':'none'})

    def _draw_net_iv(self, ht, net_iv):
        """Plot a net interval, which includes the vertex identifier label and a different color."""
        self._draw_interval(ht, net_iv, color=self.colours["blue"], mid_label=self._tr_graph.get_identifier(net_iv.nb_set))

    def gens(self, filename,start=None,**kwargs):
        self.net_intervals(self._tr_graph.net_ivs_below(start=start), filename,**kwargs)

    def graph(self, filename, edge_labels=None):
        """Draw the transition graph as specified by tr_graph (defaults to the internal transition graph).
        If edge_labels is true, also plots the edge index on the graph (it's pretty ugly).
        The vertices of the graph are colored as follows:
        - light gray: vertices that are not in any loop class or the essential class
        - red: vertices that are in the essential class, or if the transition graph is incomplete, the vertices which have out degree 0
        - other colours: the loop classes
        """
        graph = self._tr_graph.g
        comps = self._tr_graph.components() # labels for the connected components
        pos = gt.draw.sfdp_layout(graph) # graph vertex and edge positioning
        max_v = comps.a.max()

        def set_color(n):
            if n == -1:
                return self.base_light_color # light gray
            elif n == 0:
                return self.highlight_color # red
            else:
                return matplotlib.colors.to_hex(self.cmap(n/max_v))

        vtx_colors = graph.new_vertex_property("string")
        for v in graph.vertices():
            vtx_colors[v] = set_color(comps[v])

        ht = 100*int(np.sqrt(graph.num_vertices()))+130
        if edge_labels == 'index':
            gt.draw.graph_draw(graph, pos, vertex_size=14,vertex_text=graph.vertex_index,vertex_fill_color=vtx_colors,
                    edge_pen_width=2, output=filename, output_size=(ht,ht),
                    edge_text=graph.edge_index,edge_text_distance=0,edge_text_parallel=False,edge_font_size=10)
        else:
            gt.draw.graph_draw(graph, pos, vertex_size=14,vertex_text=graph.vertex_index,vertex_fill_color=vtx_colors,
                    edge_pen_width=2, output=filename, output_size=(ht,ht))

        # create legend with entries
        # - square with colour for essential class
        # - square with colour for non-loop
        # - rectangle with viridis (or any colormap) colours, for loop classes
        #  legend_elements = [matplotlib.lines.Line2D([0], [0], color=base_color, marker='o', label='Base Vertex'),
                           #  matplotlib.lines.Line2D([0], [0], color=highlight_color, marker='o', label='Essential Class')]
        #  matplotlib.pyplot.gca().legend(handles=legend_elements)


