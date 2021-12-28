import os
import pandas as pd
import math

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

    for module in ["openmpi", "intel"]:
        top_directory = "csv/" + module

        for current_directory in os.listdir(top_directory):

            current_directory_path = os.path.join(
                top_directory, current_directory)

            dataFrames_current_directory = []
            for index, filename in enumerate(os.listdir(current_directory_path)):

                filename_path = os.path.join(
                    current_directory_path, filename)

                headers, lines = removeHeadersFromCSV(filename_path)

                df = pd.read_csv(filename_path)

                insertHeadersToCSV(headers, filename_path)
                dataFrames_current_directory.append(df)

            counter = index + 1

            sum_dataFrame = dataFrames_current_directory[0]
            for df in dataFrames_current_directory[1:]:
                sum_dataFrame['t[usec]'] += df['t[usec]']
                sum_dataFrame['Mbytes/sec'] += df['Mbytes/sec']

            mean_dataFrame = sum_dataFrame

            mean_dataFrame['t[usec]'] = (
                mean_dataFrame['t[usec]'] / counter).round(2)

            mean_dataFrame['Mbytes/sec'] = (
                mean_dataFrame['Mbytes/sec'] / counter).round(2)

            dev = 0
            for i in range(0, counter):
                dev += math.sqrt((dataFrames_current_directory[i]['t[usec]'][0]
                                  - mean_dataFrame['t[usec]'][0])**2)
            std_dev = round(math.sqrt(dev / (counter - 1)), 2)

            mean_result_path = os.path.join("csv/mean/" + module, filename)
            mean_result_path = mean_result_path[:-6] + mean_result_path[-4:]

            mean_dataFrame.to_csv(mean_result_path, sep=",", index=False)

            insertHeadersToCSV(headers, mean_result_path)
            print("\nmean result saved in: " + mean_result_path,
                  "\nstandard deviation of lambda: ", std_dev)
