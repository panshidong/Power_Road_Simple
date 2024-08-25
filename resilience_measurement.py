#input is a strategy, like a repair sequence
#loop
#spend time, repair a node, assume repair is instant
#accumulate and compute resilience triangle
#re-check the functioning bus nodes if the repaired node is a bus node
#change the network parameters and re-run the model
#until sequence end

#repair sequence is some sequence of all broken links and bus
#use heuristic to find out the optimum solution
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

def load_disrupted_scenatio(broken_buses,broken_links):
    unfunctional_nodes = delete_buses(broken_buses)
    capacity_adjustment(Org_network,Network1,broken_links,0) #delete link equal to change capacity into 0
    power_to_road(unfunctional_nodes,Network1,Network2)   #This will edit the capacity of roadway link due to traffic light
    files=[]
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S-%f")[:23]
    os.makedirs(backup_dir, exist_ok=True)
    #backup_filepath=backup_dir+timestamp
    #shutil.copy2(Network2, backup_filepath)
    #files.append(timestamp)
    if os.path.exists('flows.txt'):
        os.remove('flows.txt')
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
    return len(get_functional_nodes(broken_buses))/33

def resilience_triangle(functionality,time):
    #求解若干个梯形面积之和
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
        current_node,current_move_time=repair_path_time('flows.txt',repair_seq[0],previous_node)
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

def heuristic_find_solution(initial_sequence,consider_interdependence):
    if len(initial_sequence) <= 1:
        raise ValueError("Initial sequence must contain more than one element.")

    if hasattr(creator, 'FitnessMin'):
        del creator.FitnessMin
    if hasattr(creator, 'Individual'):
        del creator.Individual
    # 创建最小化适应度类
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMin)

    toolbox = base.Toolbox()

    # 定义个体的生成规则
    toolbox.register("individual", tools.initIterate, creator.Individual, lambda: random.sample(initial_sequence, len(initial_sequence)))

    # 定义种群的生成规则
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    # 定义适应度函数
    def eval_one_max(individual):
        single_run_time_0=datetime.now()
        fitness = resilience_evaluation(individual)[0]
        single_run_time=datetime.now()-single_run_time_0
        print(f"Individual: {individual}, Fitness: {fitness}, Duration: {single_run_time}")  # 调试输出
        return (fitness,)

    toolbox.register("evaluate", eval_one_max)

    # 注册遗传算法的操作函数
    if consider_interdependence==True:
        toolbox.register("mate", cxOrdered)
        toolbox.register("mutate", mutShuffleIndexes, indpb=0.2)
    else:
        toolbox.register("mate", cxOrderedGrouped)
        toolbox.register("mutate", mutShuffleIndexesGrouped, indpb=0.2)
    toolbox.register("select", tools.selTournament, tournsize=3)

    # 初始化种群
    population = toolbox.population(n=50)
    print("Initial population:")  # 调试输出
    for ind in population[:5]:  # 只打印前5个个体
        print(ind)
    
    # 定义遗传算法的参数
    NGEN = 20  # 迭代次数
    CXPB = 0.5  # 交叉概率
    MUTPB = 0.2  # 突变概率
    
    # 记录日志
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", lambda x: sum([a[0] for a in x])/len(x))
    stats.register("min", min)
    stats.register("max", max)

    # 运行遗传算法
    algorithms.eaSimple(population, toolbox, cxpb=CXPB, mutpb=MUTPB, ngen=NGEN, 
                        stats=stats, verbose=True)
    
    # 找到最优个体
    best_ind = tools.selBest(population, 1)[0]
    with open("ans.txt", 'w') as file:
        file.writelines(str(best_ind))
    print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values))
    return best_ind

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

myind=heuristic_find_solution(sequence,True)
#myind=sequence #this is used for debug

run_end_time=datetime.now()
duration=run_end_time - run_start_time
#seperate final back up nets with others

result_opt, road_opt, power_opt, time_opt,net_files=resilience_evaluation(myind)
#for the best solution, draw the resilience triangle
plot_triangles_seperate(road_opt,power_opt,time_opt,result_folder+'opt')
plot_triangle_tot(road_opt,power_opt,time_opt,result_folder+'opt')
with open(result_folder+'output.txt', 'w') as f:
    print("This is optimal considering interdependence", file=f)
    print(myind, file=f)
    print("run duration: " + str(duration), file=f)
    print("total complement resilience(not average): ", result_opt, file=f)
    print("road resilience: ", road_opt, file=f)
    print("power resilience: ", power_opt, file=f)
    print("time steps: ", time_opt, file=f)
    print()


run_start_time=datetime.now()
myind=heuristic_find_solution(sequence,False)
run_end_time=datetime.now()
result_opt, road_opt, power_opt, time_opt,net_files=resilience_evaluation(myind)
#for the best solution, draw the resilience triangle
plot_triangles_seperate(road_opt,power_opt,time_opt,result_folder)
plot_triangle_tot(road_opt,power_opt,time_opt,result_folder)
with open(result_folder+'output.txt', 'a') as f:
    print("This is optimal NOT considering interdependence", file=f)
    print(myind, file=f)
    print("run duration: " + str(duration), file=f)
    print("total complement resilience(not average): ", result_opt, file=f)
    print("road resilience: ", road_opt, file=f)
    print("power resilience: ", power_opt, file=f)
    print("time steps: ", time_opt, file=f)

