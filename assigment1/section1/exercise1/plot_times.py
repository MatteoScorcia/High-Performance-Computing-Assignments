import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def T_comm(message_size, c, bandwidth, num_procs):
    return (c*2 + (message_size / (bandwidth * 1000000)*2)) * num_procs


if __name__ == '__main__':
    df = pd.read_csv("ring_times.csv")

    statistics = df.groupby(['#Processors'], as_index="False")[
        'Time'].agg([('average', 'mean'), ('standard_deviation', 'std')])

    plt.rcParams["figure.figsize"] = (12, 10)
    plt.scatter(statistics.index, statistics['average'],
                label="average ring time")
    plt.xticks(statistics.index)
    plt.xlabel("#Processors")
    plt.ylabel("t[sec]")
    plt.legend(loc="upper left")
    plt.savefig("plot_ring_times.jpg")

    df = pd.read_csv("runtimes.csv")

    statistics = df.groupby(['#Processors'], as_index="False")[
        'Time'].agg([('average', 'mean'), ('standard_deviation', 'std')])

    plt.scatter(statistics.index, statistics['average'],
                label="average runtime")
    plt.xticks(statistics.index)
    plt.xlabel("#Processors")
    plt.ylabel("t[sec]")

    x = np.arange(1, 25)
    c = 0.41 / 1000000
    bandwidth = 5460
    message_size = 4
    y = [T_comm(message_size, c, bandwidth, num_proc) for num_proc in x]
    plt.plot(x, y, 'k-', label="network model (#Processors <= 24)")

    x = np.arange(25, 49)
    c = 0.99 / 1000000
    bandwidth = 12060
    message_size = 4
    y = [T_comm(message_size, c, bandwidth, num_proc) for num_proc in x]
    plt.plot(x, y, 'k--', label="network model (#Processors > 24)")

    plt.legend(loc="upper left")
    plt.savefig("plot_times.jpg")
    plt.clf()

    speedup = [statistics['average'].iloc[0] /
               statistics['average'].iloc[i] for i in np.arange(0, 48)]

    efficiency = [speedup[num_proc] / (num_proc+1)
                  for num_proc in np.arange(0, 48)]

    plt.plot(np.arange(0, 48), efficiency, label="Speedup Efficiency")
    plt.xlabel("#Processors")
    plt.ylabel("Efficiency")
    plt.legend(loc="upper right")
    plt.savefig("efficiency.jpg")
