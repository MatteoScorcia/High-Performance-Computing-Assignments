import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.core.indexes.base import Index


def compute_C(L, k):
    return L**2 * k * 2 * 8


def compute_T_c(C, B, k, T_l):
    return (C/B) + (k*T_l)


def compute_MLUPs(L, N, T_s, T_c):
    return (L**3 * N) / (T_s + T_c)


if __name__ == '__main__':

    # first exercise
    results_core = []
    results_socket = []
    L = 1200
    T_s = 15.24

    for mapping in ["core", "socket"]:
        if mapping == "core":
            B = 6220 * 2 * 1000  # (full duplex, B/s)
            T_l = 0.2 / 1000000  # (s)
        else:
            B = 5460 * 2 * 1000
            T_l = 0.41 / 1000000
        for N in [4, 8, 12]:
            if (N == 4):
                k = 2
            else:
                k = 3

            C = compute_C(L, k)
            T_c = compute_T_c(C, B, k, T_l)

            res = compute_MLUPs(L, N, T_s, T_c) / 1000000
            results_core.append(
                res) if mapping == "core" else results_socket.append(res)

    plt.rcParams["figure.figsize"] = (12, 10)

    x = [4, 8, 12]
    y = [results_core[0], results_core[1], results_core[2]]

    plt.scatter(x, y, label='expected performance')
    y_measured = [451.812471379, 890.246539591, 1326.29329586]
    plt.scatter(x, y_measured, marker="^", label="measured performance")

    plt.ylabel('P(L,N) [MLUPs/sec]')
    plt.xlabel('#Processes')
    plt.legend(loc="upper left")

    plt.savefig("plots/exercise1/map-core-plot.jpg")
    plt.clf()

    diff = np.array(y_measured) - np.array(y)

    df = pd.DataFrame()
    df['#Processors'] = [4, 8, 12]
    df['expected performance [MLUPs/sec]'] = y
    df['measured performance [MLUPs/sec]'] = y_measured
    df['difference in performance [MLUPs/sec]'] = diff
    df.to_csv("csv/exercise1/map-core.csv", sep=",", index=False)

    print("exercise 1 - map by core:")
    print("expected -> ", y)
    print("computed -> ", y_measured)
    print("diff -> ", diff)
    print()

    y = [results_socket[0], results_socket[1], results_socket[2]]

    plt.scatter(x, y, label='expected performance')
    y_measured = [450.740879615, 891.285363587, 1336.25329837]
    plt.scatter(x, y_measured, marker="^", label="measured performance")

    plt.ylabel('P(L,N) [MLUPs/sec]')
    plt.xlabel('#Processes')
    plt.legend(loc="upper left")

    plt.savefig("plots/exercise1/map-socket-plot.jpg")
    plt.clf()

    diff = np.array(y_measured) - np.array(y)

    df = pd.DataFrame()
    df['#Processors'] = [4, 8, 12]
    df['expected performance [MLUPs/sec]'] = y
    df['measured performance [MLUPs/sec]'] = y_measured
    df['difference in performance [MLUPs/sec]'] = diff
    df.to_csv("csv/exercise1/map-socket.csv", sep=",", index=False)

    print("exercise 1 - map by socket")
    print("expected -> ", y)
    print("measured -> ", y_measured)
    print("diff -> ", diff)
    print()

    # second exercise
    results = []
    L = 1200
    T_s = 15.24

    B = 12060 * 2 * 1000
    T_l = 0.99 / 1000000
    k = 3

    for N in [12, 24, 48]:
        C = compute_C(L, k)
        T_c = compute_T_c(C, B, k, T_l)

        res = compute_MLUPs(L, N, T_s, T_c) / 1000000
        results.append(res)

    x = [12, 24, 48]
    y = [results[0], results[1], results[2]]

    plt.scatter(x, y, label='expected performance')
    y_measured = [1335.38757673, 2655.37926409, 5097.28061230]
    plt.scatter(x, y_measured, marker="^", label="measured performance")

    plt.ylabel('P(L,N) [MLUPs/sec]')
    plt.xlabel('#Processes')
    plt.legend(loc="upper left")

    plt.savefig("plots/exercise2_plot.jpg")
    plt.clf()

    diff = np.array(y_measured) - np.array(y)

    df = pd.DataFrame()
    df['#Processors'] = [12, 24, 48]
    df['expected performance [MLUPs/sec]'] = y
    df['measured performance [MLUPs/sec]'] = y_measured
    df['difference in performance [MLUPs/sec]'] = diff
    df.to_csv("csv/exercise2.csv", sep=",", index=False)

    print("exercise 2:")
    print("expected -> ", y)
    print("measured -> ", y_measured)
    print("diff -> ", np.array(y_measured) - np.array(y))
    print()

    results = []

    L = 1200
    T_s = 22.07

    B = 12020 * 2 * 1000
    T_l = 1.39 / 1000000
    k = 3

    for N in [12, 24, 48]:
        C = compute_C(L, k)
        T_c = compute_T_c(C, B, k, T_l)

        res = compute_MLUPs(L, N, T_s, T_c) / 1000000
        results.append(res)

    x = [12, 24, 48]
    y = [results[0], results[1], results[2]]

    plt.scatter(x, y, label='expected performance')
    y_measured = [918.746449966, 1779.42073415, 2750.75611807]
    plt.scatter(x, y_measured, marker="^", label="measured performance")

    plt.ylabel('P(L,N) [MLUPs/sec]')
    plt.xlabel('#Processes')
    plt.legend(loc="upper left")

    plt.savefig("plots/exercise3_plot.jpg")
    plt.clf()

    diff = np.array(y_measured) - np.array(y)

    df = pd.DataFrame()
    df['#Processors'] = [12, 24, 48]
    df['expected performance [MLUPs/sec]'] = y
    df['measured performance [MLUPs/sec]'] = y_measured
    df['difference in performance [MLUPs/sec]'] = diff
    df.to_csv("csv/exercise3.csv", sep=",", index=False)

    print("exercise 3:")
    print("expected -> ", y)
    print("measured -> ", y_measured)
    print("diff -> ", np.array(y_measured) - np.array(y))
    print()
