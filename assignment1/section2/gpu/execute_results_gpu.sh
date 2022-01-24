#! /bin/sh

#remember to activate your python environment with all dependecies installed before the run

python process_mean_csv_gpu.py
python fit_and_plot_gpu.py
python tabulate_results_gpu.py