from road_util import capacity_adjustment,calculate_shortest_path_cost
import os
def bus_lookup(bus):
    # this gives the location of each bus
    references = {
        1: (),
        2: (),
        3: (),
        4: (),
        5: (),
        6: (),
        7: (),
        8: (),
        9: (),
        10: (),
        11: (3,1),
        12: (),
        13: (),
        14: (),
        15: (1,3),
        16: (),
        17: (20,22),
        18: (),
        19: (),
        20: (),
        21: (),
        22: (9,10),
        23: (),
        24: (),
        25: (),
        26: (),
        27: (),
        28: (20,22),
        29: (),
        30: (),
        31: (),
        32: (22,20),
        33: ()
    }
    return references[bus]

def bus_loc(bus):
    # this gives the location of each bus
    references = {
        1: 1,
        2: 2,
        3: 3,
        4: 3,
        5: 3,
        6: 4,
        7: 4,
        8: 5,
        9: 5,
        10: 7,
        11: 7,
        12: 8,
        13: 8,
        14: 9,
        15: 10,
        16: 11,
        17: 12,
        18: 13,
        19: 14,
        20: 15,
        21: 16,
        22: 17,
        23: 18,
        24: 19,
        25: 20,
        26: 21,
        27: 22,
        28: 23,
        29: 24,
        30: 24,
        31: 20,
        32: 20,
        33: 27
    }
    return references[bus]

def repair_path_time(f, repaired,O):
    #this assumes each repair requires dispatch rather than routing
    # O is the real location on roadmap, which is determined in the previous step,  repair need to be found
    if isinstance(repaired,int):
        if bus_loc(repaired)>0:
            D=bus_loc(repaired)
            return D,calculate_shortest_path_cost('s.txt',O,D)
        else:
            return 0
    else:
        D=repaired
        front_node_cost=calculate_shortest_path_cost('s.txt',O,D[0])
        rear_node_cost=calculate_shortest_path_cost('s.txt',O,D[1])
        if front_node_cost<rear_node_cost:
            return D[0],front_node_cost
        else:
            return D[1],rear_node_cost

    """
    flow file looks like this
    (1,2) 4495.677114 6.000817 
    4495.677114 
    """
    return None

def power_to_road(unfunctional_nodes, input_file, output_file):
    adj_links=[]
    for bus in unfunctional_nodes:                                 # not an empty tuple
        if bus_lookup(bus):
            adj_links.append(bus_lookup(bus))
    capacity_adjustment(input_file, output_file, adj_links,0.5)    # traffic light reduce capacity by 50%
    return 0

