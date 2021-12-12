import pandas as pd
import os
import matplotlib.pyplot as plt
from numpy import arange
from scipy.optimize import curve_fit
import re


def objective(x_data, c, bandwidth):
    return c + (x_data / bandwidth)


if __name__ == '__main__':

    directory = "1csv"
    for filename in os.listdir(directory):
        full_path = os.path.join(directory, filename)
        print(full_path)
        df = pd.read_csv(full_path)

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
        plot_full_path = os.path.join("plots/", plot_filename)
        plt.savefig(plot_full_path)
        plt.clf()

        df['t[usec] computed'] = [round(number, 2) for number in y_estimated]
        Mbytes_comp = [round(number, 2) for number in (x_data / y_estimated)]
        df['Mbytes/sec computed'] = Mbytes_comp

        result_full_path = os.path.join("results/", filename)
        df.to_csv(result_full_path, sep="/t", index=False)

        str = os.path.splitext(filename).__str__()
        chunks = re.split('[-]', str)

        # //todo aggiungi le 3 righe di commenti
        # "#mpirun  --map-by $type --mca pml $i --mca btl self,$j -np 2 --report-bindings ./IMB-MPI1 PingPong -msglog 28"
        # list of nodes involved
        # lamba -> df['t[usec] computed'][0], bandwith -> df['Mbytes/sec computed'][-1] (computed by fitting data)
