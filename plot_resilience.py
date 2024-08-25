import matplotlib.pyplot as plt
import numpy as np

def plot_triangles_seperate(road_resilience, power_resilience,time,loc):
    # Adjusting the resilience values to start at the provided initial resilience value and end at resilience 1.0

    road_resilience_adjusted = road_resilience + [1.0]
    power_resilience_adjusted = power_resilience + [1.0]

     # Calculate cumulative time
    cumulative_time = np.cumsum(time)

    # Adjusting the cumulative time to include time 0
    cumulative_time_adjusted = [0] + list(cumulative_time)

    # Plotting Resilience Triangles with the correctly adjusted values
    plt.figure(figsize=(14, 8))

    # Road Resilience Triangle
    plt.plot(cumulative_time_adjusted, road_resilience_adjusted, '-o', label='Road Resilience', color='blue')
    plt.fill_between(cumulative_time_adjusted, 1, road_resilience_adjusted, color='blue', alpha=0.3)

    # Power Resilience Triangle
    plt.plot(cumulative_time_adjusted, power_resilience_adjusted, '-o', label='Power Resilience', color='green')
    plt.fill_between(cumulative_time_adjusted, 1, power_resilience_adjusted, color='green', alpha=0.3)

    # Labels and Title
    plt.xlabel('Cumulative Time')
    plt.ylabel('Resilience')
    plt.ylim(0,1)
    plt.title('Resilience Triangles for Road and Power Systems (Adjusted Initial Value)')
    plt.legend()

    plt.savefig(loc+'plot1.png')
    # Display the plot
    #plt.grid(True)
    #plt.show()

def plot_triangle_tot(road_resilience, power_resilience,time,loc):
    # Adjusting the resilience values to start at the provided initial resilience value and end at resilience 1.0

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
    plt.plot(cumulative_time_adjusted, total_resilience, '-o', label='Total Resilience', color='purple')
    plt.fill_between(cumulative_time_adjusted, 1, total_resilience, color='purple', alpha=0.3)

    # Labels and Title
    plt.xlabel('Cumulative Time')
    plt.ylabel('Resilience')
    plt.ylim(0,1)
    plt.title('Total Resilience of Road and Power Systems')
    plt.legend()
    plt.savefig(loc+'plot2.png')

    # Display the plot
    #plt.grid(True)
    #plt.show()

