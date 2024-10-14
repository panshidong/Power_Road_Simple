import matplotlib.pyplot as plt
import numpy as np

def plot_triangles_seperate(road_resilience, power_resilience,time,loc):
    # Adjusting the resilience values to start at the provided initial resilience value and end at resilience 1.0
    plt.rcParams.update({'font.size': 22,'font.serif':'Times New Roman',})
    road_resilience_adjusted = road_resilience + [1.0]
    power_resilience_adjusted = power_resilience + [1.0]

     # Calculate cumulative time
    cumulative_time = np.cumsum(time)

    # Adjusting the cumulative time to include time 0
    cumulative_time_adjusted = [0] + list(cumulative_time)

    # Plotting Resilience Triangles with the correctly adjusted values
    plt.figure(figsize=(14, 8))

    # Road Resilience Triangle
    plt.plot(cumulative_time_adjusted, road_resilience_adjusted, '-o', label='Highway', color='blue')
    plt.fill_between(cumulative_time_adjusted, 1, road_resilience_adjusted, color='blue', alpha=0.3)

    # Power Resilience Triangle
    plt.plot(cumulative_time_adjusted, power_resilience_adjusted, '-o', label='Power', color='green')
    plt.fill_between(cumulative_time_adjusted, 1, power_resilience_adjusted, color='green', alpha=0.3)

    # Labels and Title
    plt.xlabel('Cumulative Time')
    plt.ylabel('Performance Level')
    plt.ylim(0,1)
    plt.title('Resilience Triangles for Highway and Power Systems')
    plt.legend()

    plt.savefig(loc+'plot1.png')
    # Display the plot
    plt.grid(True)
    plt.show()

def plot_triangle_tot(road_resilience, power_resilience,time,loc):
    # Adjusting the resilience values to start at the provided initial resilience value and end at resilience 1.0
    plt.rcParams.update({'font.size': 22,'font.serif':'Times New Roman',})
    road_resilience_adjusted = road_resilience + [1.0]
    power_resilience_adjusted = power_resilience + [1.0]

     # Calculate cumulative time
    cumulative_time = np.cumsum(time)

    # Adjusting the cumulative time to include time 0
    cumulative_time_adjusted = [0] + list(cumulative_time)
    # Calculating total resilience by taking the average of road and power resilience at each time step
    total_resilience = [(r + p) / 2 for r, p in zip(road_resilience_adjusted, power_resilience_adjusted)]

    # Plotting Total Resilience
    plt.figure(figsize=(14, 8))

    # Total Resilience Curve
    plt.plot(cumulative_time_adjusted, total_resilience, '-o', label='Weighted  total', color='purple')
    plt.fill_between(cumulative_time_adjusted, 1, total_resilience, color='purple', alpha=0.3)

    # Labels and Title
    plt.xlabel('Cumulative Time')
    plt.ylabel('Performance Level')
    plt.ylim(0,1)
    plt.title('Weighted Total Resilience of Highway and Power Systems',pad=20)
    plt.legend()
    plt.savefig(loc+'plot2.png')

    # Display the plot
    plt.grid(True)
    plt.show()


def plot_triangles_compare2(resilience1, resilience2,time1,time2,desciprtion1,description2,loc):
    # Adjusting the resilience values to start at the provided initial resilience value and end at resilience 1.0
    plt.rcParams.update({'font.size': 22,'font.serif':'Times New Roman',})
    resilience1_adjusted = resilience1 + [1.0]
    resilience2_adjusted = resilience2 + [1.0]

     # Calculate cumulative time
    cumulative_time1 = np.cumsum(time1)
    cumulative_time2 = np.cumsum(time2)

    # Adjusting the cumulative time to include time 0
    cumulative_time1_adjusted = [0] + list(cumulative_time1)
    cumulative_time2_adjusted = [0] + list(cumulative_time2)

    # Plotting Resilience Triangles with the correctly adjusted values
    plt.figure(figsize=(14, 8))

    # Resilience1 Triangle 1
    plt.plot(cumulative_time1_adjusted, resilience1_adjusted, '-o', label=desciprtion1, color='blue')
    plt.fill_between(cumulative_time1_adjusted, 1, resilience1_adjusted, color='blue', alpha=0.3)

    # Resilience Triangle 2
    plt.plot(cumulative_time2_adjusted, resilience2_adjusted, '-o', label=description2, color='green')
    plt.fill_between(cumulative_time2_adjusted, 1, resilience2_adjusted, color='green', alpha=0.3)

    # Labels and Title
    plt.xlabel('Cumulative Time')
    plt.ylabel('Performance')
    plt.ylim(0,1)
    plt.title('Resilience Triangles for '+desciprtion1+ ' and '+ description2, pad=20)
    plt.legend()

    plt.savefig(loc+'plot1.png')
    # Display the plot
    plt.grid(True)
    plt.show()

def plot_triangles_compare3(resilience1, resilience2,resilience3,time1,time2,time3,desciprtion1,description2,description3,loc):
    # Adjusting the resilience values to start at the provided initial resilience value and end at resilience 1.0
    plt.rcParams.update({'font.size': 22,'font.serif':'Times New Roman',})
    resilience1_adjusted = resilience1 + [1.0]
    resilience2_adjusted = resilience2 + [1.0]
    resilience3_adjusted = resilience3 + [1.0]

     # Calculate cumulative time
    cumulative_time1 = np.cumsum(time1)
    cumulative_time2 = np.cumsum(time2)
    cumulative_time3 = np.cumsum(time3)

    # Adjusting the cumulative time to include time 0
    cumulative_time1_adjusted = [0] + list(cumulative_time1)
    cumulative_time2_adjusted = [0] + list(cumulative_time2)
    cumulative_time3_adjusted = [0] + list(cumulative_time3)

    # Plotting Resilience Triangles with the correctly adjusted values
    plt.figure(figsize=(14, 8))

    # Resilience1 Triangle 1
    plt.plot(cumulative_time1_adjusted, resilience1_adjusted, '-o', label=desciprtion1, color='red')
    plt.fill_between(cumulative_time1_adjusted, 1, resilience1_adjusted, color='red', alpha=0.3)

    # Resilience Triangle 2
    plt.plot(cumulative_time2_adjusted, resilience2_adjusted, '-o', label=description2, color='yellow')
    plt.fill_between(cumulative_time2_adjusted, 1, resilience2_adjusted, color='yellow', alpha=0.3)

    plt.plot(cumulative_time3_adjusted, resilience3_adjusted, '-o', label=description3, color='green')
    plt.fill_between(cumulative_time3_adjusted, 1, resilience3_adjusted, color='green', alpha=0.3)

    # Labels and Title
    plt.xlabel('Cumulative Time')
    plt.ylabel('Performance')
    plt.ylim(0,1)
    plt.title('Resilience Triangles for '+desciprtion1+ ' and '+ description2+ ' and '+ description3,pad=20)
    plt.legend()
    # Display the plot
    plt.grid(True)
    plt.show()
