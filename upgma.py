import sys
import graphviz
import heapq as hq

class Edge:
    """
    class Edge to represent the distance between two nodes
    (or taxa in a phylogenetic tree application)

    Fields:
        distance : float - distance (physical or evolutionary) between nodes
        node1 : str - node1 and node2 ordered st that node1 < node2 lexographically
        node2 : str
    """
  # constructor
    def __init__(self, dist, n1, n2):
        """ Constructor for class Edge to represent the distance between two nodes.

            :param dist: float distance between the two nodes
            :param n1: string rep of first node
            :param n2: string rep of second node
        """
        self.distance = dist
        self.node1 = n1
        self.node2 = n2
  
   # override the comparison operator
    def __lt__(self, nxt):
        """ Returns true if first Edge has lower distance
            If tied for distance, returns true for first lexicographically less than second
        """
        if self.distance == nxt.distance:
            return str(order_pair(self.node1, self.node2)) <\
                   str(order_pair(nxt.node1, nxt.node2))
        else: 
            return self.distance < nxt.distance

class NodeInfo:
    """
        class NodeInfo to contain data accumulated about a node/cluster
        used as values in node_dict mapping node names with commas + parens
        to NodeInfo objects
    """
    #constructor
    def __init__(self, size, level, name):
        """ Constructor for class to contain data accumulated about a node

            :param size: int - number of nodes contained within this node/cluster
            :param level: int - upgma height for DOT diagram
            :param name: str - shorthand name for the node (ie no parentheses and commas)
        """
        self.size = size
        self.level = level
        self.name = name

def parse_args():
    """
    Parses input system arguments:
        1 - filepath for distance matrix
        2 - filepath to write result to

    Parameters: NA

    Returns: 
        distances : list(str) representing pairwise distances
        write_to_file : str is filepath to write result to
    """
    distances = open(sys.argv[1],'r').readlines()
    distances = [dist for dist in distances if dist != "\n" and dist != " "]
    write_to_file = sys.argv[2]
    return distances, write_to_file

def order_pair(node1, node2):
    """
    Consumes 2 nodes, returns them as a tuple in lexicographic order

        :param node1: - str
        :param node2: - str

    Returns
        ordered_pair : tuple(str, str)
    """
    if node1 < node2:
        ordered_pair = (node1, node2)
    else:
        ordered_pair = (node2, node1)
    return ordered_pair

def instantiate_distances(distances):
    """ Consumes list of distances, returns data structures to be used in neighbor-joining.

        :param distances: list(str) representing pairwise distances

    Returns: 
        distance_heap : heapq of Edges sorted by distance
        node_dict : dict[str, NodeInfo] mapping each cluster to number of member nodes
        distance_dict : dict[(str, str), float] mapping tuple of nodes to distance
                                                between the members of the pair

    """
    distance_heap, node_dict, dist_dict = [], {}, {}

    for dist in distances:
        nodes_dist = dist.split()
        node1, node2, distance = nodes_dist[0], nodes_dist[1], float(nodes_dist[2])
        if node1 < node2: 
            new_edge = Edge(distance, node1, node2)
        else: 
            new_edge = Edge(distance, node2, node1)
        dist_dict[order_pair(node1, node2)] = distance
        ONE_NODE_START = 1
        ZERO_UPGMA_HEIGHT = 0
        node_dict[node1] = NodeInfo(ONE_NODE_START, ZERO_UPGMA_HEIGHT, node1)
        node_dict[node2] = NodeInfo(ONE_NODE_START, ZERO_UPGMA_HEIGHT, node2)
        hq.heappush(distance_heap, new_edge)
    return distance_heap, node_dict, dist_dict

def calc_distance(node, union, node_dict, dist_dict):
    """
    Calculates the distance between a given node and the newly unioned node
    
    Parameters: 
        :param node: - str rep of any node
        :param union: - Edge that has just been added to the graph
        :param node_dict: - dict[node: NodeInfo] mapping string node to it's NodeInfo
        :param dist_dict: - dict[(node1, node2): distance] mapping pair to the distance between them

    Returns:
        distance : float
    """
    node_to_ui = dist_dict[order_pair(union.node1, node)]
    node_to_uj = dist_dict[order_pair(union.node2, node)]
    num_items_i, num_items_j = node_dict[union.node1].size, node_dict[union.node2].size
    coeff_i = num_items_i / (num_items_i + num_items_j)
    coeff_j = num_items_j / (num_items_i + num_items_j)
    distance = (coeff_i * node_to_ui) + (coeff_j * (node_to_uj)) 
    return distance

def update_heap_and_dict(cluster, dist_heap, node_dict, dist_dict):
    """
    Updates given distance_heap to reflect union of nodes in given cluster

    Parameters: 
        cluster : Edge containing nodes being clustered
        dist_heap : heap(Edge)
        node_dict : dict[node: int]
        dist_dict : dict[(node1, node2): float]

    Returns: 
        updated distance_heap and distance_dict

    Side affects: 
        Adds distances to the given distances_heap between union of closest_nodes
        and every non-member node in the nodes_set
        Removes elements which contain individual members of new union
    """
    valid_edge = lambda dist :\
            dist.node1 != cluster.node1\
        and dist.node1 != cluster.node2\
        and dist.node2 != cluster.node1\
        and dist.node2 != cluster.node2
    #update dist_heap by taking out edges to now-unioned nodes
    new_dist_heap = []
    [hq.heappush(new_dist_heap, dist) for dist in list(dist_heap) if valid_edge(dist)]

    # add new distances from cluster to heap and dict
    nodes = [node for node in node_dict.keys() if node != cluster.node1 and node != cluster.node2]
    for node in nodes:
        new_dist = calc_distance(node, cluster, node_dict, dist_dict)
        cluster_name = "(" + cluster.node1 + "," + cluster.node2 + ")"
        if node < cluster_name:
            new_edge = Edge(new_dist, node, cluster_name)
        else:
            new_edge = Edge(new_dist, cluster_name, node)
        hq.heappush(new_dist_heap, new_edge)
        dist_dict[order_pair(node, cluster_name)] = new_dist

    return new_dist_heap, dist_dict

def update_node_dict(cluster, node_dict):
    """
    Updates given node_dict to reflect union of given cluster

    Parameters: 
        cluster : Edge containing nodes being clustered
        node_dict : dict[str: NodeInfo(size, level)]

    Returns: 
        updated node_dict
    """
    # save components of new dict entry
    new_union = "(" + cluster.node1 + "," + cluster.node2 + ")"
    new_node_count = node_dict[cluster.node1].size + node_dict[cluster.node2].size
    new_node_level = max(node_dict[cluster.node1].level, node_dict[cluster.node2].level) + 1
    new_node_name = node_dict[cluster.node1].name + node_dict[cluster.node2].name
    # add entry for new unioned node
    node_dict[new_union] = NodeInfo(new_node_count, new_node_level, new_node_name)

    return node_dict

def add_graph_cluster(dot_output, cluster, node_dict):
    """
    Adds edges to given Graph object to represent given cluster

    Parameters: 
        dot_output : Graph - graph object to add edges to
        cluster : Edge - cluster being formed
        node_dict : dict(str, NodeInfo) - maps names of current nodes to NodeInfo
    
    Returns: 
        updated Graph object with edges from latest cluster
    """
    node1, node2 = cluster.node1, cluster.node2
    node1_level, node2_level  = node_dict[node1].level, node_dict[node2].level

    cluster_name = "(" + cluster.node1 + "," + cluster.node2 + ")"
    union_level = node_dict[cluster_name].level
    union_name = node_dict[cluster_name].name + str(union_level)

    dot_output.edge(node_dict[node1].name + str(node1_level), union_name)
    dot_output.edge(node_dict[node2].name + str(node2_level), union_name)

    return dot_output

def main(distances, output_file):
    """
    Outputs the cluster pattern in terminal and writes graphviz graph to given dot file

    Parameters: 
        distances : list(str) representing all original graph edges
        output_file : str - filepath to write output to

    Returns: NA

    Side affects: 
    Prints clustering pattern to terminal
    Writes a phylogenetic treeg in DOT format to input .dot file 
    """
    distance_heap, node_dict, distance_dict = instantiate_distances(distances)
    dot_output = graphviz.Graph('tree')
    while len(node_dict) > 1:
        cluster = hq.heappop(distance_heap)
        distance_heap, distance_dict = update_heap_and_dict(cluster, distance_heap, node_dict, distance_dict)
        node_dict = update_node_dict(cluster, node_dict)
        dot_output = add_graph_cluster(dot_output, cluster, node_dict)
        # delete entries for old nodes - used to break while-loop
        del node_dict[cluster.node1]
        del node_dict[cluster.node2]

    [sys.stdout.write(key) for key in node_dict.keys()]
    dot_output.save(output_file)

if __name__ == "__main__":
    distances, write_to_file = parse_args()
    main(distances, write_to_file)