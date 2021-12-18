#! /bin/sh

conda activate ml #my virtual environment

python process_mean_csv.py
python fit_and_plot.py
python tabulate_results.py