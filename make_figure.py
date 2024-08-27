import random
import networkx as nx
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
