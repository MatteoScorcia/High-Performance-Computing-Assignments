import os
import pandas as pd

# take mean of csv and make 1 result in same folder to be processed by fit_and_plot.py


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
    print()
