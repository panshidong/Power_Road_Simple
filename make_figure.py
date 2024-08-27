import matplotlib.pyplot as plt
import networkx as nx
import random

# Define the tree-like network structure
connections = {
    1: [2],
    2: [3, 19],
    3: [4, 23],
    4: [5],
    5: [6],
    6: [7, 26],
    7: [8],
    8: [9, 21],
    9: [10],
    10: [11],
    11: [12],
    12: [13, 22],
    13: [14],
    14: [15],
    15: [16],
    16: [17],
    17: [18],
    18: [33],
    19: [20],
    20: [21],
    21: [],
    22: [],
    23: [24],
    24: [25],
    25: [29],
    26: [27],
    27: [28],
    28: [29],
    29: [30],
    30: [31],
    31: [32],
    32: [],
    33: []
}

# Create a directed graph
G = nx.DiGraph()
for node, successors in connections.items():
    for successor in successors:
        G.add_edge(node, successor)

# Define a custom function to create a hierarchical tree layout
def hierarchy_pos(G, root=None, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5):
    pos = _hierarchy_pos(G, root, width, vert_gap, vert_loc, xcenter)
    return pos

def _hierarchy_pos(G, root, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5, pos=None, parent=None, parsed=[]):
    if pos is None:
        pos = {root: (xcenter, vert_loc)}
    else:
        pos[root] = (xcenter, vert_loc)
        
    children = list(G.neighbors(root))
    if not isinstance(G, nx.DiGraph) and parent is not None:
        children.remove(parent)  
            
    if len(children) != 0:
        dx = width / len(children) 
        nextx = xcenter - width/2 - dx/2
        for child in children:
            nextx += dx
            pos = _hierarchy_pos(G, child, width=dx, vert_gap=vert_gap, 
                                 vert_loc=vert_loc-vert_gap, xcenter=nextx,
                                 pos=pos, parent=root, parsed=parsed)
    return pos

# Use the custom function to get the layout
pos = hierarchy_pos(G, root=1)

# Randomly select three nodes to color red
random.seed(42)
red_nodes = random.sample(list(G.nodes), 3)

# Define node colors: red for selected nodes, blue for others
node_colors = ['red' if node in red_nodes else 'cornflowerblue' for node in G.nodes]

# Draw the nodes and edges with the specified nodes in red
plt.figure(figsize=(14, 12))
nx.draw(G, pos, 
        with_labels=True, 
        node_color=node_colors, 
        node_size=800, 
        font_size=12, 
        font_color='white',
        font_weight='bold', 
        edge_color='darkgray', 
        linewidths=1.5, 
        arrows=True, 
        arrowstyle='-|>', 
        arrowsize=20)

# Add a title with a larger font size
plt.title("Enhanced Hierarchical Tree-like Network Structure with Highlighted Nodes", fontsize=16, fontweight='bold')

# Display the final figure
plt.show()

# Print the red nodes
print("Red nodes are:", red_nodes)
