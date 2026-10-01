import random
def get_functional_nodes(broken_nodes):
    # Define the topology of the IEEE 33-bus system   this function is broken, use the below one
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
    # FIX (TaskC_OD branch, 2026-09-03): a bus is functional iff neither it nor any
    # upstream bus is broken. `connections` maps parent -> children on the radial
    # feeder, so a broken bus takes down itself and every descendant. The previous
    # implementation recursed over *children* with a shared `visited` cache and
    # therefore did not propagate outages downstream (see
    # backup_prefix_powerbug_20260903/BACKUP_RECORD.md for the analysis).
    broken = set(int(b) for b in broken_nodes)
    unfunctional = set(broken)
    stack = list(broken)
    while stack:
        current = stack.pop()
        for child in connections.get(current, []):
            if child not in unfunctional:
                unfunctional.add(child)
                stack.append(child)

    functional_nodes = set(connections.keys()) - unfunctional
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
