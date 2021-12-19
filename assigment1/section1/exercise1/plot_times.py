import pandas as pd
import os
import matplotlib.pyplot as plt
from numpy import arange
from scipy.optimize import curve_fit

if __name__ == '__main__':
    df = pd.read_csv("ring_times.csv")

    plt.rcParams["figure.figsize"] = (12, 10)
    plt.scatter(df['#Processors'], df['MeanTime'] *
                1000000, label="time vs #processors")
    plt.xticks(df['#Processors'])
    plt.xlabel("#Processors")
    plt.ylabel("t[usec]")
    plt.savefig("plot_times.jpg")
    plt.clf()
