import pandas as pd
import os
import matplotlib.pyplot as plt
from tabulate import tabulate
import csv


if __name__ == '__main__':

    full_path = "table.csv"

    rows = []
    with open(full_path, "r") as file:
        csvreader = csv.reader(file)
        titles = next(csvreader)
        for row in csvreader:
            rows.append(row)

        res = tabulate(rows, titles, numalign="right")
        res_full_path = os.path.join(full_path[:-3]+"txt")

        with open(res_full_path, "w+") as file:
            file.writelines(res)
