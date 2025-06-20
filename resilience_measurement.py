'''
input is a strategy, like a repair sequence
loop
spend time, repair a node, assume **a constant repair time**
accumulate and compute resilience triangle
re-check the functioning bus nodes if the repaired node is a bus node
change the network parameters and re-run the model
until sequence end
repair sequence is some sequence of all broken links and bus
use heuristic to find out the optimum solution ***for***:
***resilience only
***equity only
***both
'''
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
import numpy as np

power_road_factor=0.5
broken_link_factor=0
cached_values={} #this will explode if number of broken items too big, may want to use some kind of rolling techniques

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
def evalequity(seq,times,criteria):
    '''
    #ASSUMPTION----------equity of power and road are given same weight
    '''
    #REVIEW HERE
    def gini_coefficient(x):
        """
        计算数组 x 的 Gini 系数。
        实现参考：StackOverflow 上的向量化实现
        :contentReference[oaicite:0]{index=0}
        """
        x = np.sort(np.array(x, dtype=float))
        n = x.size
        if n == 0:
            return 0.0
        index = np.arange(1, n+1)
        return (np.sum((2*index - n - 1) * x)) / (n * np.sum(x))

    # 1. 累计时间
    accum_time = list(itertools.accumulate(times))  # :contentReference[oaicite:1]{index=1}

    # 2. 按类型分组收集停运时间
    bus_times = []
    road_times = []
    for item, t in zip(seq, accum_time):             # :contentReference[oaicite:2]{index=2}
        if isinstance(item, int):                    # :contentReference[oaicite:3]{index=3}
            bus_times.append(t)
        else:
            road_times.append(t)
    #3. 根据 criteria 选择不同的衡量函数
    if criteria == 'var':
        # 方差：数据离散程度的经典度量 :contentReference[oaicite:4]{index=4}
        power_eq = np.var(bus_times) if bus_times else 0
        road_eq  = np.var(road_times) if road_times else 0

    elif criteria == 'mad':
        # 平均绝对偏差（Mean Absolute Deviation）:contentReference[oaicite:5]{index=5}
        def mad(x):
            x = np.array(x, dtype=float)
            return np.mean(np.abs(x - np.mean(x))) if x.size>0 else 0
        power_eq = mad(bus_times)
        road_eq  = mad(road_times)

    elif criteria == 'cv':
        # 变异系数（Coeff. of Variation）= std/mean :contentReference[oaicite:6]{index=6}
        def cv(x):
            x = np.array(x, dtype=float)
            return (np.std(x) / np.mean(x)) if x.size>0 and np.mean(x)!=0 else 0
        power_eq = cv(bus_times)
        road_eq  = cv(road_times)

    elif criteria == 'gini':
        # Gini 系数 :contentReference[oaicite:7]{index=7}
        power_eq = gini_coefficient(bus_times)
        road_eq  = gini_coefficient(road_times)
    return (power_eq+road_eq)/2


def resilience_evaluation(repair_seq):
    #This gives resilience "triangle"
    repair_seq=repair_seq.copy()
    resilience_road=[]
    resilience_power=[]
    time=[]
    net_file_names=[]
    previous_node=13             #this is the assumed initial depot
    repair_time=0     
    while len(repair_seq)>0:
        broken_buses=[]
        broken_links=[]
        for item in repair_seq:
            if isinstance(item,int):
                broken_buses.append(item)
                repair_time=20            #assumed bus bar repair time
            else:
                broken_links.append(item)
                repair_time=10           # assumed link clearance time
        #cache process
        brokens=tuple(repair_seq)
        if brokens in cached_values:
            (current_resilience_road,
            current_resilience_power,
            total_time,
            current_node) = cached_values[brokens]
        else:
            #take road output and read travel time, give number
            current_resilience_road=eval_road_resilience(broken_buses,broken_links)
            #take power output and give number
            current_resilience_power=eval_power_resilience(broken_buses)
            #repair and continue
            current_node,current_move_time=repair_path_time('s.txt',repair_seq[0],previous_node)
            total_time = current_move_time + repair_time
            #store in to cached
            cached_values[brokens] = (current_resilience_road, current_resilience_power,total_time, current_node)
        resilience_road.append(current_resilience_road)
        resilience_power.append(current_resilience_power)
        time.append(total_time)         
        broken_buses=[bus for bus in broken_buses if bus!=repair_seq[0]]
        broken_links=[link for link in broken_links if link!=repair_seq[0]]
        #set up for next loop
        previous_node=current_node
        repair_seq.pop(0)
        #net_file_names.append(load_disrupted_scenatio(broken_buses,broken_links))
    full_resilience = resilience_triangle(resilience_road,time)+resilience_triangle(resilience_power,time)
    criteria_list = ['var', 'mad', 'cv', 'gini'] 
    equity = {c: evalequity(repair_seq, time, c) for c in criteria_list}
    return full_resilience, resilience_road,resilience_power,time,net_file_names,equity

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

def find_solution_all(initial_sequence,focus):
    result_opt=-1
    print(resilience_evaluation([17, (9, 10), (11, 14), 15, 11, 32, 28]))
    for seq in itertools.permutations(initial_sequence):
        seq=list(seq)
        #REVIEW change the criteria here
        if focus=='Resilience':
            result=resilience_evaluation(seq)[0]
        else:
            result=resilience_evaluation(seq)[5][focus]
        if result_opt > resilience_result or result_opt<0:
            result_opt=resilience_result
            best_seq=seq
    return best_seq


def heuristic_find_solution(initial_sequence,consider_interdependence):
    start_time=time.time()
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
        #FIXME------------------------------------------
        #if use quity then [0] change to [5]
        fitness = resilience_evaluation(individual)[0]
        #fitness = resilience_evaluation(individual)[5]['gini']
        single_run_time=datetime.now()-single_run_time_0
        #print(f"Individual: {individual}, Fitness: {fitness}, Duration: {single_run_time}")  # 调试输出
        return (fitness,)

    toolbox.register("evaluate", eval_one_max)

    # 注册遗传算法的操作函数
    if consider_interdependence==True:
        toolbox.register("mate", cxOrdered)
        toolbox.register("mutate", mutShuffleIndexes, indpb=0.2)
    else:
        toolbox.register("mate", cxOrderedGrouped)
        toolbox.register("mutate", mutShuffleIndexesGrouped, indpb=0.2)
    toolbox.register("select", tools.selTournament, tournsize=5)

    # 初始化种群
    population = toolbox.population(n=100)
    print("Initial population:")  # 调试输出
    for ind in population[:5]:  # 只打印前5个个体
        print(ind)
    
    # 定义遗传算法的参数
    NGEN = 50  # 迭代次数
    CXPB = 0.5  # 交叉概率
    MUTPB = 0.2  # 突变概率
    
    # 记录日志
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", lambda x: sum([a[0] for a in x])/len(x))
    stats.register("min", min)
    stats.register("max", max)
    print("definitions:"+str(time.time()-start_time))
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
#REVIEW change the initial sequence to avoid comparison
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

def run_model(sequence,bool_stream,result_folder,message,Scenario,plot_control,focus):
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
        #REVIEW----------------this is the solution method
        #myind=heuristic_find_solution(sequence,bool_stream)
        myind=find_solution_all(sequence,focus)
    #myind=sequence #this is used for debug

    run_end_time=datetime.now()
    duration=run_end_time - run_start_time
    #seperate final back up nets with others

    result_opt, road_opt, power_opt, time_opt,net_files,equity_results=resilience_evaluation(myind)
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
        if focus=="Resilience":
            print("Equity:Variance ", evalequity(sequence,time_opt,'var'),file=f)
            print("Equity:mad ", evalequity(sequence,time_opt,'mad'),file=f)
            print("Equity:cv ", evalequity(sequence,time_opt,'cv'),file=f)
            print("Equity:gini ", evalequity(sequence,time_opt,'gini'),file=f)
        else:
            print("Equity ", evalequity(sequence,time_opt,focus),file=f)
        print()

    return myind

#IDEA the following senario are tested, default is randon sequence, then there is the resilience maximum without interdependency,
#IDEA  then there is resilience maximum with interdependencies,showing interdependency is important but
#IDEA then there is maximum equity with different criteria, with interdependency, comparing with resilience maximum showing equity is also important
#IDEA then there is maximum equity with no interdependency showing interdependencies are also important when it comes to equity considerations
# IDEA  then there are also some less important sensitivities that i dont intend to do, listed as is like the previous study 
#sequence=[11,17,15,(9,10),28,32,(11,14)]
s_resilience="Resilience"
run_model(sequence,True,result_folder,"This is random",'evalrand',True,s_resilience)   #consider default sequence as random
run_model(sequence,True,result_folder,"This is optimal considering interdependence",'opt',True,s_resilience) #optimal resilience considering interdependency
run_model(sequence,True,result_folder,"This is optimal considering interdependence",'opt',True,'var') #optimal resilience considering interdependency
run_model(sequence,True,result_folder,"This is optimal considering interdependence",'opt',True,'mad') #optimal resilience considering interdependency
run_model(sequence,True,result_folder,"This is optimal considering interdependence",'opt',True,'cv') #optimal resilience considering interdependency
run_model(sequence,True,result_folder,"This is optimal considering interdependence",'opt',True,'gini') #optimal resilience considering interdependency


powers_only=[11,17,15,28,32]
roads_only=[(9,10),(11,14)]
equity_measures=['var','mad','cv','gini']
for r in equity_measures:
    run_model(sequence,True,result_folder,"This is optimal considering interdependence",'opt',True,r) #optimal resilience considering interdependency
    power_ans=run_model(powers_only,True,result_folder,"This is optimal Power only",'optPower',False,r)
    roads_ans=run_model(roads_only,True,result_folder,"This is optimal Road only",'optRoad',False,r)
    roads_ans=run_model(roads_ans+power_ans,True,result_folder,"This is Road priority",'evalRoadPriority',True,r)
    roads_ans=run_model(power_ans+roads_ans,True,result_folder,"This is Power priority",'evalPowerPriority',True,r)

run_model(sequence,False,result_folder,"This is optimal NOT considering interdependence",'',False)


"""
'''
Sensitivity design
Sensitivity #1: 
increase the number of broken links/nets
'''
SENS1_sequence=[11,17,15,(9,10),28,32,(11,14),(15,22),(2,6),6,24]
run_model(SENS1_sequence,True,result_folder,"This is random of SENSITIVITY #1",'evalrandSENS1',True)
run_model(SENS1_sequence,True,result_folder,"This is SENSITIVITY #1",'SENS1',True)

powers_only=[11,17,15,28,32,6,24]
roads_only=[(9,10),(11,14),(15,22),(2,6)]
power_ans=run_model(powers_only,True,result_folder,"This is optimal Power only SENS1",'optPowerSENS1',False)
roads_ans=run_model(roads_only,True,result_folder,"This is optimal Road only SENS1",'optRoadSENS1',False)
roads_ans=run_model(roads_ans+power_ans,True,result_folder,"This is Road priority SENS1",'evalRoadPrioritySENS1',True)
roads_ans=run_model(power_ans+roads_ans,True,result_folder,"This is Power priority SENS1",'evalPowerPrioritySENS1',True)

'''
Sensitivity#2:
move the connection points around

'''
run_model(sequence,True,result_folder,"This is SENSITIVITY #3",'SENS2',True)
powers_only=[11,17,15,28,32]
roads_only=[(9,10),(11,14)]
power_ans=run_model(powers_only,True,result_folder,"This is optimal Power only SENS2",'optPowerSENS2',False)
roads_ans=run_model(roads_only,True,result_folder,"This is optimal Road only SENS2",'optRoadSENS2',False)
roads_ans=run_model(roads_ans+power_ans,True,result_folder,"This is Road priority SENS2",'evalRoadPrioritySENS2',True)
roads_ans=run_model(power_ans+roads_ans,True,result_folder,"This is Power priority SENS2",'evalPowerPrioritySENS2',True)

'''
Sensitivity #3
different harm level for broken net or power fail (interdependency level)

'''
power_road_factor=0.3 #the lower of this the more severe the damage is
run_model(sequence,True,result_folder,"This is SENSITIVITY #3",'SENS3',True)
powers_only=[11,17,15,28,32]
roads_only=[(9,10),(11,14)]
power_ans=run_model(powers_only,True,result_folder,"This is optimal Power only SENS3",'optPowerSENS3',False)
roads_ans=run_model(roads_only,True,result_folder,"This is optimal Road only SENS3",'optRoadSENS3',False)
roads_ans=run_model(roads_ans+power_ans,True,result_folder,"This is Road priority SENS3",'evalRoadPrioritySENS3',True)
roads_ans=run_model(power_ans+roads_ans,True,result_folder,"This is Power priority SENS3",'evalPowerPrioritySENS3',True)
power_road_factor=0.5

'''
Sensitivity #4
Interdependency Pattern like where the interdependenct location is

'''
run_model(sequence,True,result_folder,"This is SENSITIVITY #4",'SENS4',True)
powers_only=[11,17,15,28,32]
roads_only=[(9,10),(11,14)]
power_ans=run_model(powers_only,True,result_folder,"This is optimal Power only SENS4",'optPowerSENS4',False)
roads_ans=run_model(roads_only,True,result_folder,"This is optimal Road only SENS4",'optRoadSENS4',False)
roads_ans=run_model(roads_ans+power_ans,True,result_folder,"This is Road priority SENS4",'evalRoadPrioritySENS4',True)
roads_ans=run_model(power_ans+roads_ans,True,result_folder,"This is Power priority SENS4",'evalPowerPrioritySENS4',True)


'''
Sensitivity #5
One directional interdependency

'''
power_road_factor=1.0
run_model(sequence,True,result_folder,"This is SENSITIVITY #5",'SENS4',True)
powers_only=[11,17,15,28,32]
roads_only=[(9,10),(11,14)]
power_ans=run_model(powers_only,True,result_folder,"This is optimal Power only SENS5",'optPowerSENS5',False)
roads_ans=run_model(roads_only,True,result_folder,"This is optimal Road only SENS5",'optRoadSENS5',False)
roads_ans=run_model(roads_ans+power_ans,True,result_folder,"This is Road priority SENS5",'evalRoadPrioritySENS5',True)
roads_ans=run_model(power_ans+roads_ans,True,result_folder,"This is Power priority SENS5",'evalPowerPrioritySENS5',True)
power_road_factor=0.5


'''
Sensitivity #X not used
Demand pattern, like 50% demand in at time 0 and gradual recover? hard to design


'''
"""