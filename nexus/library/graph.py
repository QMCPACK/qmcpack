##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  graph.py                                                          #
#    Interface with pydot to make images of directed graphs.         #
#                                                                    #
#  Content summary:                                                  #
#    Graph                                                           #
#      Wrapper class for pydot functionality                         #
#                                                                    #
#====================================================================#


from generic import obj
from developer import DevBase,unavailable
try:
    from pydot import Dot,Node,Edge,Cluster as pdCluster
except ImportError:
    Dot,Node,Edge,pdCluster = unavailable('Dot','Node','Edge','pdCluster')
#end try


# websites
#  http://www.graphviz.org/Home.php
#  http://pythonhaven.wordpress.com/tag/pydot/
#  http://jseabold.net/blog/2012/02/making-graphical-models-with-pydot.html


class Graph(DevBase):

    pydot_type = Dot

    def __init__(self,*args,**kwargs):
        # some args are
        #   graph_name = 'G'
        #   graph_type in ('graph','digraph')
        self.graph = self.pydot_type(*args,**kwargs)
        self.name  = self.graph.get_name()
        self.nodes = obj()
        self.edges = obj()
        self.subgraphs = obj()
    #end def __init__


    def set_node_defaults(self,*args,**kwargs):
        #  shape, height, width, fontsize
        self.graph.set_node_defaults(*args,**kwargs)
    #end def set_node_defaults


    def get_node(self,name):
        return self.nodes[name]
    #end def get_node


    def add_node(self,*args,**kwargs):
        # some args are
        #   name = 'N'
        #   texlbl = 'N_{n,m}'
        if len(args)==1 and isinstance(args[0],Node):
            node = args[0]
        else:
            node = Node(*args,**kwargs)
        #end if
        name = node.get_name()
        self.nodes[name] = node
        self.graph.add_node(node)
    #end def add_node


    def add_nodes(self,*nodes):
        if len(nodes)==1 and isinstance(nodes[0],list):
            nodes = nodes[0]
        #end if
        for node in nodes:
            if isinstance(node,Node):
                self.add_node(node)
            else:
                self.add_node(**node)
            #end if
        #end for
    #end def add_nodes


    def del_node(self,name):
        del self.nodes[name]
        self.graph.del_node(name)
    #end def del_node


    def del_nodes(self,*names):
        for name in names:
            self.del_node(name)
        #end for
    #end def del_nodes


    def get_edge(self,name):
        return self.edges[name]
    #end def get_edge


    def add_edge(self,*args,**kwargs):
        if len(args)==1 and isinstance(args[0],Edge):
            edge = args[0]
        else:
            edge = Edge(*args,**kwargs)
        #end if
        src = edge.get_source()
        dst = edge.get_destination()
        self.edges[src,dst] = edge
        self.graph.add_edge(edge)
    #end def add_edge


    def add_edges(self,*edges):
        if len(edges)==1 and isinstance(edges[0],list):
            edges = edges[0]
        #end if
        for edge in edges:
            if isinstance(edge,Edge):
                self.add_edge(edge)
            else:
                self.add_edge(**edge)
            #end if
        #end for
    #end def add_edges


    def del_edge(self,name):
        del self.edges[name]
        self.graph.del_edge(name)
    #end def del_edge


    def del_edges(self,*names):
        for name in names:
            self.del_edge(name)
        #end for
    #end def del_edges


    def incorporate(self,graph):
        if not isinstance(graph,Graph):
            self.error('graph must be of type Graph\n  you provided: '+graph.__class__.__name__)
        #end if
        self.add_nodes(*graph.nodes)
        self.add_edges(*graph.edges)
    #end def incorporate


    def get_subgraph(self,name):
        return self.subgraphs[name]
    #end def get_subgraph

    def add_subgraph(self,graph):
        #self.incorporate(graph)
        self.subgraphs[graph.name] = graph
        self.graph.add_subgraph(graph)
    #end def add_subgraph


    def add_subgraphs(self,*graphs):
        if len(graphs)==1 and isinstance(graphs[0],list):
            graphs = graphs[0]
        #end if
        for graph in graphs:
            self.add_subgraph(graph)
        #end for
    #end def add_subgraphs


    def del_subgraph(self,name):
        sg = self.get_subgraph(name)
        self.del_edges(*sg.edges.list())
        self.del_nodes(*sg.nodes.list())
        del self.subgraphs[name]
    #end def del_subgraph


    def del_subgraphs(self,*names):
        for name in names:
            self.del_subgraph(name)
        #end for
    #end def del_subgraphs


    def write_png(self,*args,**kwargs):
        self.graph.write_png(*args,**kwargs)
    #end def write_png


    def write_dot(self,*args,**kwargs):
        self.graph.write_dot(*args,**kwargs)
    #end def write_dot


    def write_pdf(self,*args,**kwargs):
        self.graph.write_pdf(*args,**kwargs)
    #end def write_pdf


    def write_svg(self,*args,**kwargs):
        self.graph.write_svg(*args,**kwargs)
    #end def write_svg


    def write_xdot(self,*args,**kwargs):
        self.graph.write_xdot(*args,**kwargs)
    #end def write_xdot
#end class Graph




class Cluster(Graph):
    pydot_type = pdCluster
#end class Cluster
