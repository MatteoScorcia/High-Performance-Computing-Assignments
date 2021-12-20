import pandas as pd
import os
import matplotlib.pyplot as plt
from numpy import arange
from scipy.optimize import curve_fit

if __name__ == '__main__':
    df = pd.read_csv("ring_times.csv")

    statistics = df.groupby(['#Processors'], as_index="False")[
        'Time'].agg([('average', 'mean'), ('standard_deviation', 'std')])
    print(statistics)
    print()

    plt.rcParams["figure.figsize"] = (12, 10)
    plt.scatter(statistics.index, statistics['average'],
                label="average time vs #processors")
    plt.xticks(statistics.index)
    plt.xlabel("#Processors")
    plt.ylabel("t[sec]")
    plt.savefig("plot_times.jpg")
    plt.clf()
