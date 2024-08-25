import matplotlib.pyplot as plt
import numpy as np
def plot_triangles_seperate(resilience1, resilience2,time1,time2):
    # Adjusting the resilience values to start at the provided initial resilience value and end at resilience 1.0

    resilience1 = resilience1 + [1.0]
    resilience2 = resilience2 + [1.0]

     # Calculate cumulative time
    cumulative_time1 = np.cumsum(time1)
    cumulative_time2 = np.cumsum(time2)

    # Adjusting the cumulative time to include time 0
    cumulative_time_adjusted1 = [0] + list(cumulative_time1)
    cumulative_time_adjusted2 = [0] + list(cumulative_time2)

    # Plotting Resilience Triangles with the correctly adjusted values
    plt.figure(figsize=(14, 8))

    # Road Resilience Triangle
    plt.plot(cumulative_time_adjusted1, resilience1,'-o', label='Resilience1', color='blue')
    plt.fill_between(cumulative_time_adjusted1, 1, resilience1, color='blue', alpha=0.3)

    # Power Resilience Triangle
    plt.plot(cumulative_time_adjusted2, resilience2, '-o', label='Resilience2', color='green')
    plt.fill_between(cumulative_time_adjusted2, 1, resilience2, color='green', alpha=0.3)

    # Labels and Title
    plt.xlabel('Cumulative Time')
    plt.ylabel('Resilience')
    plt.ylim(0,1)
    plt.title('Resilience Triangles for Road and Power Systems (Adjusted Initial Value)')
    plt.legend()
    plt.grid(True)
    plt.show()

road_resilience1 =  [0.4737320321335463, 0.47371965441681574, 0.4736659980722235, 0.4737261827951427, 0.4938443259966816, 0.5066794705606961, 0.6309284976026788]
power_resilience1 =  [0.696969696969697, 0.7575757575757576, 0.8181818181818182, 0.8787878787878788, 0.9393939393939394, 1.0, 1.0]
time1 =[3.117423, 82.889342, 81.235567, 66.105362, 77.86277100000001, 31.819896, 20.174809]
resilience1=[]
for i in range(len(road_resilience1)):
    resilience1.append((road_resilience1[i]+power_resilience1[i])/2)


road_resilience2=  [0.4737731975168442, 0.7220992593835087, 0.8643246746142932, 0.9324386137591868, 0.9332979866107854, 0.9498025898906521, 0.9741399291762326]
power_resilience2=  [0.696969696969697, 0.696969696969697, 0.696969696969697, 0.7575757575757576, 0.8181818181818182, 0.8787878787878788, 0.9393939393939394]
time2 =  [21.484966, 21.849974000000003, 48.630033999999995, 29.495787, 21.448567, 3.031529, 38.856482]
resilience2=[]
for i in range(len(road_resilience2)):
    resilience2.append((road_resilience2[i]+power_resilience2[i])/2)

plot_triangles_seperate(resilience1,resilience2,time1,time2)