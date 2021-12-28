import pandas as pd
import os
import matplotlib.pyplot as plt
from numpy import arange
from scipy.optimize import curve_fit
import re
from tabulate import tabulate


def objective(x_data, c, bandwidth):
    return c + (x_data / bandwidth)


def removeHeadersFromCSV(full_path):
    with open(full_path, "r") as f:
        lines = f.readlines()

    with open(full_path, "w") as f:
        f.writelines(lines[3:])

    return lines[0:3], lines[3:]


def insertHeadersToCSV(headers, result_full_path):
    with open(result_full_path, "r") as f:
        lines = f.readlines()

    with open(result_full_path, "w") as f:
        for index, head in enumerate(headers):
            lines.insert(index, head)
        f.writelines(lines)


if __name__ == '__main__':

    for module in ["openmpi", "intel"]:
        directory = "csv/mean/" + module

        for filename in os.listdir(directory):
            full_path = os.path.join(directory, filename)

            headers, lines = removeHeadersFromCSV(
                full_path)

            df = pd.read_csv(full_path)

            insertHeadersToCSV(headers, full_path)

            x_data, y = df['#bytes'].values, df['t[usec]'].values
            tick_list = list(range(len(x_data - 1)))
            popt, _ = curve_fit(objective, x_data, y)
            a, b = popt
            plt.rcParams["figure.figsize"] = (12, 10)

            plt.scatter(x_data, y, label='measured data')
            y_estimated = objective(x_data, a, b)
            plt.plot(x_data, y_estimated, '--', color='red',
                     label='least squares fitting')

            c = y[0]
            bandwidth = df['Mbytes/sec'][-1:]

            T_comm = [c + (element / bandwidth) for element in x_data]
            plt.plot(x_data, T_comm, '--', color='blue',
                     label='simple communication model')
            plt.ylabel('t[usecon]')
            plt.xlabel('#bytes')
            plt.xscale('log', base=2)
            plt.yscale('log', base=10)
            plt.legend(loc="upper left")

            plot_filename = os.path.splitext(filename)[0]+'.jpg'
            plot_full_path = os.path.join("plots/" + module, plot_filename)
            plt.savefig(plot_full_path)
            plt.clf()

            df['t[usec] (computed)'] = [round(number, 2)
                                        for number in y_estimated]
            Mbytes_comp = [round(number, 2)
                           for number in (x_data / y_estimated)]
            df['Mbytes/sec (computed)'] = Mbytes_comp

            result_full_path = os.path.join("results/" + module, filename)
            df.to_csv(result_full_path, sep=",", index=False)

            headers[2] = '#lambda[usec] (computed) -> ' + \
                str(df['t[usec] (computed)'].iloc[0]) + \
                ', bandwidth[Mbytes/sec] (computed) -> ' + \
                str(df['Mbytes/sec (computed)'].iloc[-1]) + '\n'

            insertHeadersToCSV(headers, result_full_path)
