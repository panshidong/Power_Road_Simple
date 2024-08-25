import random
def get_functional_nodes(broken_nodes):
    # Define the topology of the IEEE 33-bus system
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
    # Function to check if a node is functional
    def is_functional(node, broken_nodes, connections, visited):
        if node in broken_nodes:
            return False
        if node in visited:
            return True
        visited.add(node)
        for neighbor in connections[node]:
            if not is_functional(neighbor, broken_nodes, connections, visited):
                return False
        return True

    # Initialize the set of functional nodes
    functional_nodes = set()
    visited = set()

    # Check each node
    for node in connections.keys():
        if is_functional(node, broken_nodes, connections, visited):
            functional_nodes.add(node)

    return functional_nodes

def delete_buses(broken_nodes):
    #all_buses = list(range(1, 34))  # Buses are numbered from 1 to 33
    #broken_nodes = random.sample(all_buses, num_buses_to_delete)
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
    def propagate_failure(broken, connections):
        affected = set(broken)
        to_check = set(broken)
        while to_check:
            current = to_check.pop()
            for neighbor in connections.get(current, []):
                if neighbor not in affected:
                    affected.add(neighbor)
                    to_check.add(neighbor)
        return affected

    # Propagate failure to connected nodes
    unfunctional_nodes = propagate_failure(broken_nodes, connections)
    functional_nodes = get_functional_nodes(unfunctional_nodes)
    
    all_nodes = set(connections.keys())
    unfunctional_nodes = unfunctional_nodes
    functional_nodes = all_nodes - unfunctional_nodes - set(broken_nodes)

    #return broken_nodes, unfunctional_nodes, functional_nodes
    return unfunctional_nodes
