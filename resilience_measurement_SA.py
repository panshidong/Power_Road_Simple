#input is a strategy, like a repair sequence
#loop
#spend time, repair a node, assume repair is instant
#accumulate and compute resilience triangle
#re-check the functioning bus nodes if the repaired node is a bus node
#change the network parameters and re-run the model
#until sequence end

#repair sequence is some sequence of all broken links and bus
#use heuristic to find out the optimum solution


#this is a test comment
from power_util import delete_buses
from power_util import get_functional_nodes
from road_util import capacity_adjustment
from road_util import eval_tot_OD_travel_time
from interdependency import power_to_road
from run_tapb import run_tapb
from interdependency import repair_path_time
from plot_resilience import plot_triangles_seperate,plot_triangle_tot
import random
from deap import base, creator, tools, algorithms
import itertools
import os
from datetime import datetime
import shutil
import time

power_road_factor=0.5
broken_link_factor=0


def load_disrupted_scenatio(broken_buses,broken_links):
    unfunctional_nodes = delete_buses(broken_buses)
    capacity_adjustment(Org_network,Network1,broken_links,broken_link_factor) #delete link equal to change capacity into 0
    power_to_road(unfunctional_nodes,Network1,Network2,power_road_factor)   #This will edit the capacity of roadway link due to traffic light
    files=[]
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S-%f")[:23]
    os.makedirs(backup_dir, exist_ok=True)
    #backup_filepath=backup_dir+timestamp
    #shutil.copy2(Network2, backup_filepath)
    #files.append(timestamp)
    if os.path.exists('s.txt'):
        os.remove('s.txt')
    run_tapb(Network2,'tap-b/net/SiouxFalls_trips.txt') # 1. fix newest folder issue, #2 this is ugly now, change to parameters input later
    #shutil.copy2('flows.txt', backup_filepath+'flows.txt')
    if os.path.exists(Network2):
        os.remove(Network2)
    return files

def eval_road_resilience(broken_buses,broken_links):
    load_disrupted_scenatio(broken_buses,broken_links)
    #total travel time full Sioux Falls network:7,475,338
    return 7475338/eval_tot_OD_travel_time() 

def eval_power_resilience(broken_buses):
    # This represents "unsatisfied demand" 
    return 1-len(delete_buses(broken_buses))/33

def resilience_triangle(functionality,time):
    #求解若干个梯形面积之和
    #implement financial measures for different weights
    complement=0
    functionality_for_triangle=functionality+[1]
    for i in range(len(functionality)):
        complement+=(1-functionality_for_triangle[i]+1-functionality_for_triangle[i+1])*time[i]/2
    return complement

def resilience_evaluation(repair_seq):
    #This gives resilience "triangle"
    repair_seq=repair_seq.copy()
    resilience_road=[]
    resilience_power=[]
    time=[]
    net_file_names=[]
    previous_node=13
    while len(repair_seq)>0:
        broken_buses=[]
        broken_links=[]
        for item in repair_seq:
            if isinstance(item,int):
                broken_buses.append(item)
            else:
                broken_links.append(item)
        #take road output and read travel time, give number
        resilience_road.append(eval_road_resilience(broken_buses,broken_links))
        #take power output and give number
        resilience_power.append(eval_power_resilience(broken_buses))
        #repair and continue
        current_node,current_move_time=repair_path_time('s.txt',repair_seq[0],previous_node)
        time.append(current_move_time) 
        broken_buses=[bus for bus in broken_buses if bus!=repair_seq[0]]
        broken_links=[link for link in broken_links if link!=repair_seq[0]]
        #set up for next loop
        previous_node=current_node
        repair_seq.pop(0)
        #net_file_names.append(load_disrupted_scenatio(broken_buses,broken_links))
    full_resilience = resilience_triangle(resilience_road,time)+resilience_triangle(resilience_power,time)
    return full_resilience, resilience_road,resilience_power,time,net_file_names

###########################################################################################
#This is for the comparison between optimal considering interdependency and repair by type
##########################################################################################
def cxOrderedGrouped(ind1, ind2):
    """执行有序交叉 (Order Crossover, OX)，确保不产生重复元素且保持元组和整数的分组顺序。"""
    # Split individuals into tuple and int groups
    tuples_ind1 = [x for x in ind1 if isinstance(x, tuple)]
    ints_ind1 = [x for x in ind1 if isinstance(x, int)]
    tuples_ind2 = [x for x in ind2 if isinstance(x, tuple)]
    ints_ind2 = [x for x in ind2 if isinstance(x, int)]
    
    # Apply order crossover to tuples and ints separately
    def order_crossover(part1, part2):
        size = len(part1)
        a, b = sorted(random.sample(range(size), 2))
        
        child1 = [None] * size
        child2 = [None] * size
        
        # Copy the crossover slice from the first parent to the first child
        child1[a:b + 1] = part1[a:b + 1]
        child2[a:b + 1] = part2[a:b + 1]
        
        # Fill the remaining positions with the other parent's elements
        fill_pos1, fill_pos2 = (b + 1) % size, (b + 1) % size
        for i in range(size):
            pos = (b + 1 + i) % size
            if part2[pos] not in child1:
                child1[fill_pos1] = part2[pos]
                fill_pos1 = (fill_pos1 + 1) % size
            if part1[pos] not in child2:
                child2[fill_pos2] = part1[pos]
                fill_pos2 = (fill_pos2 + 1) % size
        
        return child1, child2
    
    # Perform order crossover for both tuples and integers
    child1_tuples, child2_tuples = order_crossover(tuples_ind1, tuples_ind2)
    child1_ints, child2_ints = order_crossover(ints_ind1, ints_ind2)
    
    # Combine tuples and ints back together
    child1 = child1_tuples + child1_ints
    child2 = child2_tuples + child2_ints
    
    return creator.Individual(child1), creator.Individual(child2)

def mutShuffleIndexesGrouped(individual, indpb):
    """执行突变操作，确保不产生重复元素且保持元组和整数的分组顺序。"""
    # Split individual into tuple and int groups
    tuples_part = [x for x in individual if isinstance(x, tuple)]
    ints_part = [x for x in individual if isinstance(x, int)]
    
    # Shuffle tuples and ints separately
    def shuffle_part(part):
        size = len(part)
        for i in range(size):
            if random.random() < indpb:
                swap_indx = random.randint(0, size - 1)
                part[i], part[swap_indx] = part[swap_indx], part[i]
        return part
    
    shuffled_tuples = shuffle_part(tuples_part)
    shuffled_ints = shuffle_part(ints_part)
    
    # Combine shuffled tuples and ints back together
    shuffled_individual = shuffled_tuples + shuffled_ints
    
    return creator.Individual(shuffled_individual),
###########################################################################################
##########################################################################################
##########################################################################################

def cxOrdered(ind1, ind2):
    """执行有序交叉 (Order Crossover, OX)，确保不产生重复元素"""
    size = len(ind1)
    a, b = sorted(random.sample(range(size), 2))
    
    child1 = [None]*size
    child2 = [None]*size
    
    # Copy the crossover slice from the first parent to the first child
    child1[a:b + 1] = ind1[a:b + 1]
    child2[a:b + 1] = ind2[a:b + 1]

    # Fill the remaining positions with the other parent's elements
    fill_pos1, fill_pos2 = (b + 1) % size, (b + 1) % size
    for i in range(size):
        pos = (b + 1 + i) % size
        if ind2[pos] not in child1:
            child1[fill_pos1] = ind2[pos]
            fill_pos1 = (fill_pos1 + 1) % size
        if ind1[pos] not in child2:
            child2[fill_pos2] = ind1[pos]
            fill_pos2 = (fill_pos2 + 1) % size

    return creator.Individual(child1), creator.Individual(child2)

def mutShuffleIndexes(individual, indpb):
    """执行突变操作，确保不产生重复元素"""
    size = len(individual)
    for i in range(size):
        if random.random() < indpb:
            swap_indx = random.randint(0, size - 1)
            individual[i], individual[swap_indx] = individual[swap_indx], individual[i]
    return creator.Individual(individual),

def heuristic_find_solution(initial_sequence, consider_interdependence, *, 
                            pop_size=30, ngen=10, cxpb=0.5, mutpb=0.2, elites=1):
    start_time = time.time()
    if len(initial_sequence) <= 1:
        raise ValueError("Initial sequence must contain more than one element.")

    # Re-create DEAP classes cleanly
    if hasattr(creator, 'FitnessMin'):
        del creator.FitnessMin
    if hasattr(creator, 'Individual'):
        del creator.Individual
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMin)

    toolbox = base.Toolbox()

    # --- Initialization: permutation of the given initial_sequence ---
    toolbox.register(
        "individual",
        tools.initIterate,
        creator.Individual,
        lambda: random.sample(initial_sequence, len(initial_sequence))
    )
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    # --- Evaluation ---
    def eval_one(individual):
        single_run_time_0 = datetime.now()
        fitness = resilience_evaluation(individual)[0]
        _single_run_time = datetime.now() - single_run_time_0
        return (fitness,)
    toolbox.register("evaluate", eval_one)

    # --- Operators ---
    if consider_interdependence is True:
        # (plain permutation)
        toolbox.register("mate", cxOrdered)                          # or tools.cxPartialyMatched
        toolbox.register("mutate", mutShuffleIndexes, indpb=0.2)
    else:
        # (group-aware permutation)
        toolbox.register("mate", cxOrderedGrouped)
        toolbox.register("mutate", mutShuffleIndexesGrouped, indpb=0.2)

    # Selection should be registered in BOTH cases
    toolbox.register("select", tools.selTournament, tournsize=2)

    # --- Create initial population ---
    population = toolbox.population(n=pop_size)
    print("Initial population:")
    for ind in population[:5]:
        print(ind)

    # --- Stats / logbook ---
    stats = tools.Statistics(lambda ind: ind.fitness.values[0])
    stats.register("avg", lambda xs: sum(xs) / len(xs))
    stats.register("min", min)
    stats.register("max", max)

    print("definitions:" + str(time.time() - start_time))

    # === Elitist eaSimple (guarantees we never lose the current best) ===
    from deap.algorithms import varAnd

    # Evaluate initial population
    invalid = [ind for ind in population if not ind.fitness.valid]
    for ind, fit in zip(invalid, map(toolbox.evaluate, invalid)):
        ind.fitness.values = fit

    # Track best-so-far (monotone)
    best_so_far_vals = []
    current_best = min(ind.fitness.values[0] for ind in population)

    # Optional: Hall of Fame if you want to hold the single best individual
    hof = tools.HallOfFame(maxsize=1)
    hof.update(population)

    # Log generation 0
    gen = 0
    fits = [ind.fitness.values[0] for ind in population]
    log_min, log_avg, log_max = min(fits), sum(fits) / len(fits), max(fits)
    current_best = min(current_best, log_min)
    best_so_far_vals.append(current_best)
    print(f"{gen}\t{len(invalid)}\t{log_avg:.6f}\t({log_min},)\t({log_max},)")

    # Main loop with elitism
    for gen in range(1, ngen + 1):
        # --- keep elites from current pop ---
        elites_list = tools.selBest(population, k=elites)

        # Variation
        offspring = varAnd(population, toolbox, cxpb, mutpb)

        # Evaluate offspring
        invalid = [ind for ind in offspring if not ind.fitness.valid]
        for ind, fit in zip(invalid, map(toolbox.evaluate, invalid)):
            ind.fitness.values = fit

        # Select next generation (make room for elites)
        population = toolbox.select(population + offspring, k=len(population) - elites)

        # Insert elites unchanged (elitism)
        population.extend(map(toolbox.clone, elites_list))

        # Update HoF and stats
        hof.update(population)
        fits = [ind.fitness.values[0] for ind in population]
        log_min, log_avg, log_max = min(fits), sum(fits) / len(fits), max(fits)

        # Best-so-far is monotone by construction now
        current_best = min(current_best, log_min)
        best_so_far_vals.append(current_best)

        # Match your printed log header layout
        print(f"{gen}\t{len(invalid)}\t{log_avg:.6f}\t({log_min},)\t({log_max},)")

    # --- Result ---
    best_ind = hof[0] if len(hof) else tools.selBest(population, 1)[0]
    with open("ans.txt", "w") as f:
        f.writelines(str(best_ind))
    print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values))

    # If you want the best-so-far series for plotting, you can return it too
    return best_ind  # or: return best_ind, best_so_far_vals

import math, random, time
from typing import List, Tuple
def simulated_annealing(
    initial_sequence: List[int],
    consider_interdependence: bool = True,
    max_iters: int = 10_000,
    T0: float = 1.0,             # initial temperature
    alpha: float = 0.995,        # geometric cooling
    Tmin: float = 1e-4,          # stop when T < Tmin
    patience: int = 2000,        # stop if no improvement for this many steps
    seed: int = None,
    verbose: bool = True,
):
    if seed is not None:
        random.seed(seed)

    current = list(initial_sequence) if not isinstance(initial_sequence, list) else initial_sequence[:]

    # --- helpers ---
    def objective(seq: List[int]) -> float:
        return resilience_evaluation(seq)[0]  # minimize

    def neighbor(seq: List[int]) -> List[int]:
        n = len(seq)
        i, j = random.sample(range(n), 2)
        new_seq = seq[:]
        new_seq[i], new_seq[j] = new_seq[j], new_seq[i]
        if not consider_interdependence and random.random() < 0.2 and n >= 3:
            k = random.randrange(n)
            a = [seq[(k + d) % n] for d in range(3)]
            a = [a[2], a[0], a[1]]
            for d in range(3):
                new_seq[(k + d) % n] = a[d]
        return new_seq

    evals = 0
    cache = {}

    def cached_obj(seq: List[int]) -> float:
        nonlocal evals
        key = tuple(seq)  # if elements can be unhashable, use repr(seq)
        if key in cache:
            return cache[key]
        val = objective(seq)
        cache[key] = val
        evals += 1
        return val

    # --- init ---
    start_time = time.time()
    current_val = cached_obj(current)
    best, best_val = current[:], current_val
    T = T0
    it = 0
    no_improve = 0

    # rolling avg of accepted current states
    window = []
    window.append(current_val)

    # running mean of all candidates tried (Welford-style, numerically stable)
    cand_count = 0
    cand_mean = 0.0

    if verbose:
        print(f"[SA] start: f={current_val:.6g}, T0={T0}, alpha={alpha}, Tmin={Tmin}, max_iters={max_iters}")

    # --- main loop ---
    while it < max_iters and T > Tmin and no_improve < patience:
        it += 1
        cand = neighbor(current)
        cand_val = cached_obj(cand)

        # update candidate running mean (accepted + rejected)
        cand_count += 1
        cand_mean += (cand_val - cand_mean) / cand_count

        delta = cand_val - current_val  # we minimize

        # accept rule
        if delta <= 0 or random.random() < math.exp(-delta / T):
            current, current_val = cand, cand_val
            window.append(current_val)  # update rolling avg on accept
            if cand_val < best_val:
                best, best_val = cand[:], cand_val
                no_improve = 0
            else:
                no_improve += 1
        else:
            no_improve += 1

        # cool
        T *= alpha

        # light progress print
        if verbose and (it % 1000 == 0 or no_improve == 0):
            elapsed = time.time() - start_time
            avg_recent = sum(window) / len(window)
            print(
                f"[SA] it={it} T={T:.4g} cur={current_val:.6g} best={best_val:.6g} "
                f"avg_recent(500)={avg_recent:.6g} cand_avg={cand_mean:.6g} "
                f"no_improve={no_improve} cache={len(cache)} elapsed={elapsed:.1f}s"
            )
            print(f"[SA] total objective evaluations: {evals}")

    if verbose:
        print(f"[SA] done at it={it}, best={best_val:.6g}, time={time.time()-start_time:.2f}s")
        print(best, best_val)

    return best

run_start_time=datetime.now()
#To be replaced by relative references
Exp_folder='Experiment/'
Org_network=Exp_folder+"SiouxFalls_net.txt"
Network1 = Exp_folder + "SiouxFalls_net_link_delete.txt"
Network2 = Exp_folder+"SiouxFalls_net_use.txt"
backup_dir=Exp_folder+"Backup_nets/"
result_folder=Exp_folder+ datetime.now().strftime("%Y-%m-%d_%H-%M-%S")+'/'
os.makedirs(result_folder, exist_ok=True)
#To be replaced by random generated ones
#broken_bus_init=[11,17]
#broken_links_init=[(8,9),(9,8),(24,21),(21,24)]
sequence=[11,17,15,(9,10),28,32,(11,14)]
#print(resilience_evaluation([9,8,6,1,6,3,3]))

"""
#####################debug session###################################
myind=[(9, 10), 28, 11, 17, 15, 32, (11, 14)]
result_opt, road_opt, power_opt, time_opt,net_files=resilience_evaluation(myind)
plot_triangles_seperate(road_opt,power_opt,time_opt,result_folder+'test')
plot_triangle_tot(road_opt,power_opt,time_opt,result_folder+'test')
with open(result_folder+'output_test.txt', 'w') as f:
    print("This is optimal considering interdependence", file=f)
    print(myind, file=f)
    print("total complement resilience(not average): ", result_opt, file=f)
    print("road resilience: ", road_opt, file=f)
    print("power resilience: ", power_opt, file=f)
    print("time steps: ", time_opt, file=f)
    print("net files: ", net_files, file=f)
    print()
exit()
"""

def run_model(sequence,bool_stream,result_folder,message,Scenario,plot_control):
    if os.path.exists('bus_location.json'):
        os.remove('bus_location.json')
    if os.path.exists('bus_to_link.json'):
        os.remove('bus_to_link.json')   
    if Scenario=='SENS4':
        shutil.copy2('SENS4_bus_to_link.json', 'bus_to_link.json')
    else:
        shutil.copy2('original_bus_to_link.json', 'bus_to_link.json')
    if Scenario=='SENS2':
        shutil.copy2('SENS4_bus_location.json', 'bus_location.json')
    else:
        shutil.copy2('original_bus_location.json', 'bus_location.json')
    run_start_time=datetime.now()
    if Scenario[:4]=='eval':
        myind=sequence
    else:
        myind=simulated_annealing(sequence,bool_stream)
        myind=heuristic_find_solution(sequence,bool_stream)
    #myind=sequence #this is used for debug

    run_end_time=datetime.now()
    duration=run_end_time - run_start_time
    #seperate final back up nets with others

    result_opt, road_opt, power_opt, time_opt,net_files=resilience_evaluation(myind)
    #for the best solution, draw the resilience triangle
    if plot_control==True:
        #for the best solution, draw the resilience triangle
        plot_triangles_seperate(road_opt,power_opt,time_opt,result_folder+Scenario)
        plot_triangle_tot(road_opt,power_opt,time_opt,result_folder+Scenario)
    with open(result_folder+'output.txt', 'a') as f:
        print(message, file=f)
        print(myind, file=f)
        print("run duration: " + str(duration), file=f)
        print("total complement resilience(not average): ", result_opt, file=f)
        print("road resilience: ", road_opt, file=f)
        print("power resilience: ", power_opt, file=f)
        print("time steps: ", time_opt, file=f)
        print("-------------------------------------------------------------------------",file=f)
        print()

    return myind

for i in range(10):
    run_model(sequence,True,result_folder,"This is optimal considering interdependence",'opt',True)


