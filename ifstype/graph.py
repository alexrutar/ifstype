import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
import numpy as np
import random

class TransitionGraph(nx.DiGraph):
    def __init__(self, nb_mgr):
        super().__init__()
        self.nb_mgr = nb_mgr

    def add(self, nb_set, ch_nb_sets):
        for ch_nb in ch_nb_sets:
            self.add_edge(nb_set, ch_nb)

    def set_labels(self, eval_func):
        self.labels = {nd:eval_func(nd) for nd in self}

    def draw(self,filename="graph.pdf"):
        fig = plt.figure()
        pos = graphviz_layout(self)
        comps = list(nx.strongly_connected_components(self))
        colour_list = plt.cm.Paired(np.linspace(0,1,len(comps)))

        for idx,component in enumerate(comps):
            c = [colour_list[idx]] * len(component)
            nx.draw_networkx_nodes(component,
                    pos,
                    with_labels=False,
                    node_color=c)
        nx.draw_networkx_edges(self,pos,with_labels=False)
        nx.draw_networkx_labels(self, pos, self.labels)
        plt.axis('off')
        plt.savefig(filename)
