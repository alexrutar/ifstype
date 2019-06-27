import matplotlib
from collections import Counter
import graph_tool as gt

class TransitionGraph:
    # take any nb_mgr in which the numbering is not changing
    def __init__(self, nb_mgr):
        # graph initialization
        self.g = gt.Graph()
        for _ in range(nb_mgr.num_nb_sets()):
            self.g.add_vertex()

        for nb_set in nb_mgr:
            transition = nb_mgr.transitions[nb_set]
            for ch_nb in transition.nb_sets():
                self.g.add_edge(nb_mgr.nb_set_type(nb_set), nb_mgr.nb_set_type(ch_nb))


    def draw(self,filename="graph.pdf"):
        #  fig = plt.figure()
        #  pos = graphviz_layout(self)
        #  comps = list(nx.strongly_connected_components(self))
        #  colour_list = plt.cm.Paired(np.linspace(0,1,len(comps)))

        #  for idx,component in enumerate(comps):
            #  c = [colour_list[idx]] * len(component)
            #  nx.draw_networkx_nodes(component,
                    #  pos,
                    #  with_labels=False,
                    #  node_color=c)
        #  nx.draw_networkx_edges(self,pos,with_labels=False)
        #  nx.draw_networkx_labels(self, pos, self.labels)
        #  plt.axis('off')
        #  plt.savefig(filename)

        pos = gt.draw.sfdp_layout(self.g)
        comp, _ = gt.topology.label_components(self.g)
        is_self_loop = gt.stats.label_self_loops(self.g,mark_only=True)

        comp_size_dict = Counter(comp[v] for v in self.g.vertices()) # dictionary which takes a component and returns the number of vertices in the same component

        def non_loop_isolated(v):
            # check if a vertex is isolated and not a loop
            if comp_size_dict[comp[v]] == 1:
                for e in v.out_edges():
                    if is_self_loop[e]:
                        return False
                return True
            else:
                return False

        for v in self.g.vertices():
            if non_loop_isolated(v):
                comp[v] = -1

        print(comp[0])

        gt.draw.graph_draw(self.g, pos, vertex_size=14,vertex_text=self.g.vertex_index,vertex_fill_color=comp, edge_pen_width=2,
                   vcmap=matplotlib.cm.gist_ncar, output=filename)
