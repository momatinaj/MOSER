import networkx as nx
import random


def save_graph(G, output_path, form='QX'):
    mapping = {}
    last = 0 if form == 'ESCAPE' else 1
    edges = []
    with open(output_path, "w") as f:
        header = f"{G.number_of_nodes()}\n" if form != 'ESCAPE' else f"{G.number_of_nodes()} {G.number_of_edges()}\n"
        f.write(header)
        if form == "QX":
            edges = list(G.edges())
            edges = sorted(edges)
            result = "\n".join([f"{edge[0]} {edge[1]}" for edge in edges])
            f.write(result)
            return
        for edge in list(G.edges()):
            src, dst = edge
            if src not in mapping:
                mapping[src] = last
                last += 1
            if dst not in mapping:
                mapping[dst] = last
                last += 1
            edges.append((mapping[src], mapping[dst]))
        edges = sorted(edges)
        result = "\n".join([f"{edge[0]} {edge[1]}" for edge in edges])
        f.write(result)
    return mapping, edges


class RandomGraph:
    def __load_graph_from_escape_format(self):
        with open(self.path_to_graph, "r") as graph_file:
            graph_file.readline()
            return nx.read_edgelist(graph_file, create_using=(nx.DiGraph if self.directed else nx.Graph))

    def find_switch_candidates(self, repeat_limit=10):  # returns two edges that can be switched with each other
        m = self.G.number_of_edges()
        edge_list = list(self.G.edges)
        # draw two random edges from G
        s1, t1 = edge_list[random.randint(0, m - 1)]
        s2, t2 = None, None
        # same node edge switch avoidance
        while repeat_limit > 0:
            s2, t2 = edge_list[random.randint(0, m - 1)]
            node_list = [s1, s2, t1, t2]
            if (s2, t1) not in edge_list and (s1, t2) not in edge_list and len(list(set(node_list))) == 4:
                if self.directed:
                    break
                elif (t1, s2) not in edge_list and (t2, s1) not in edge_list:
                    break
            s2 = None
            repeat_limit -= 1
        candidate_1 = [s1, t1]
        candidate_2 = [s2, t2]
        return candidate_1, candidate_2

    def __init__(self, path_to_graph, directed=False):
        self.path_to_graph = path_to_graph
        self.directed = directed
        self.G = self.__load_graph_from_escape_format()

    def __single_random_switch(self):  # find two candidate edges and perform a switch
        e1, e2 = self.find_switch_candidates()
        e3 = [e1[0], e2[1]]
        e4 = [e2[0], e1[1]]
        if None in e1 or None in e2:
            print("unsuccessful switch!")
            return False
        self.G.remove_edges_from([e1, e2])
        self.G.add_edges_from([e3, e4])
        return True

    def mix_graph(self, number_of_switches=0):  # perform a mixing based on random switches
        if number_of_switches == 0:
            nos = self.G.number_of_edges() * 5
        else:
            nos = number_of_switches
        i = 0
        while i < nos:
            try:
                if self.__single_random_switch():
                    i += 1
                else:
                    print("trying again....")

            except Exception as e:
                print(e)

    def get_area_around_the_edge(self, edge,
                                 depth,
                                 nodelistRet=False):
        from collections import deque
        d = deque()
        visited = {edge[0]: True, edge[1]: True}
        neighbourhood = set(list(self.G.neighbors(edge[0])) + list(self.G.neighbors(edge[1])))
        for node in neighbourhood:
            d.append((node, 2))
            visited[node] = True
        while len(d) > 0:
            front = d.popleft()
            u = front[0]
            level = front[1]
            if level > depth:
                break
            for v in self.G.neighbors(u):
                if not visited.get(v, False):
                    d.append((v, level + 1))
                    visited[v] = True
        node_list = list(set(visited.keys()))

        if not nodelistRet:
            return self.G.subgraph(nodes=node_list)
        else:
            return node_list
