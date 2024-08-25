import random
import heapq
def capacity_adjustment(input_file, output_file,links,adj_factor):
    with open(input_file, 'r') as file:
        lines = file.readlines()
    
    # Initialize output lines with the metadata
    output_lines = lines[:8]  # Assuming metadata ends at line 8
    links_start = 8  # The line where links data start

    for line in lines[links_start:]:
        if line.strip().startswith('~') or line.strip().startswith(';'):
            continue  # Skip headers and comments

        parts = line.strip().split()
        if len(parts) < 10:
            continue  # Skip invalid lines

        init_node = int(parts[0])
        term_node = int(parts[1])
        capacity = float(parts[2])    # read link data from each line

        #if this line need to be modified, modify it
        for link in links:
            if adj_factor<0.1:
                if ((init_node == link[0] and term_node == link[1]) or (init_node == link[1] and term_node == link[0])):
                    parts[4]=f"{9999}"
            elif ((init_node == link[0] and term_node == link[1]) or (init_node == link[1] and term_node == link[0])):
                capacity = capacity*adj_factor

        parts[2] = f"{capacity:.8f}"
        output_lines.append('\t'.join(parts) + ' \n')

    with open(output_file, 'w') as file:
        file.writelines(output_lines)

def eval_tot_OD_travel_time():
    import csv
    
    f = "flows.txt"

    tstt = 0.0
    with open(f, 'r') as file:
        for line in file:
            # 去掉行首尾的空白字符，并检查是否为空行
            line = line.strip()
            if not line:
                continue  # 跳过空行
            
            # 检查是否是 (o,d) 行
            if line.startswith('('):
                # 将行拆分为三个部分
                parts = line.split()
                if len(parts) != 3:
                    continue  # 跳过不符合格式的行
                
                # 提取 tstt 值
                tstt += float(parts[1])*float(parts[2])
                

    return tstt

def calculate_shortest_path_cost(file_path, start, end):
    # 读取流量和花费数据
    def read_flows_and_costs(file_path):
        flows_and_costs = {}
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                if line.strip():
                    parts = line.split()
                    if len(parts) == 3:
                        node1, node2 = map(int, parts[0][1:-1].split(','))
                        flow = float(parts[1])
                        cost = float(parts[2])
                        flows_and_costs[(node1, node2)] = cost
        return flows_and_costs

    # Dijkstra算法计算最短路径
    def dijkstra(graph, start, end):
        pq = [(0, start)]  # 优先队列存储（当前总花费，当前节点）
        visited = set()  # 已访问节点集合
        min_cost = {start: 0}  # 起始节点到各节点的最小花费

        while pq:
            current_cost, current_node = heapq.heappop(pq)

            if current_node in visited:
                continue

            visited.add(current_node)

            if current_node == end:
                return current_cost

            for (u, v), cost in graph.items():
                if u == current_node and v not in visited:
                    next_cost = current_cost + cost
                    if v not in min_cost or next_cost < min_cost[v]:
                        min_cost[v] = next_cost
                        heapq.heappush(pq, (next_cost, v))

        return float('inf')  # 如果无路径，返回无穷大

    # 读取文件并计算最短路径花费
    flows_and_costs = read_flows_and_costs(file_path)
    shortest_path_cost = dijkstra(flows_and_costs, start, end)
    return shortest_path_cost

