import pandas as pd
import os
import matplotlib.pyplot as plt
from numpy import arange
from scipy.optimize import curve_fit

if __name__ == '__main__':
    df = pd.read_csv("ring_times_2.csv")

    plt.rcParams["figure.figsize"] = (12, 10)
    plt.scatter(df['#Processors'], df['Time'], label="time vs #processors")
    plt.xticks(df['#Processors'])
    plt.xlabel("#Processors")
    plt.ylabel("t[sec]")
    plt.savefig("plot_times.jpg")
    plt.clf()
