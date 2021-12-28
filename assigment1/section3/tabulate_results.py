import os
from tabulate import tabulate
import csv


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

    for subfolder in ["exercise1", ""]:
        directory = "csv/" + subfolder

        for filename in os.listdir(directory):
            full_path = os.path.join(directory, filename)

            if not os.path.isdir(full_path):

                rows = []
                with open(full_path, "r") as file:
                    csvreader = csv.reader(file)
                    titles = next(csvreader)
                    for row in csvreader:
                        rows.append(row)

                    res = tabulate(rows, titles, numalign="right")
                    tabular_full_path = os.path.join(
                        "tabular/" + subfolder, filename[:-3]+"txt")

                    with open(tabular_full_path, "w+") as file:
                        file.writelines(res)
