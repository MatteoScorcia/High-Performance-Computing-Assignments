import matplotlib.pyplot as plt
import pandas as pd

if __name__ == '__main__':
  df = pd.read_csv("strong_scalability.csv")
  processes = df['#mpi processes'].values
  runtimes = df['total runtime'].values
  speedup = [runtimes[0] / x for x in runtimes]

  print(speedup)

  plt.scatter(processes, speedup, label='strong scaling')
  plt.plot(processes, speedup)
  plt.plot(processes, processes, label='linear speedup')
  plt.ylabel('speedup')
  plt.xlabel('#mpi processes')
  plt.legend(loc="upper left")
  plt.savefig("strong_scaling_plot.jpeg")
  plt.clf()

  efficiency = [speedup[i]/processes[i] for i in range(0,len(processes))]
  plt.scatter(processes, efficiency, label='strong efficiency')
  plt.plot(processes, efficiency)
  plt.plot(processes, [1 for _ in processes], label='ideal strong efficiency')
  plt.ylabel('efficiency')
  plt.xlabel('#mpi processes')
  plt.legend(loc="lower left")
  plt.savefig("strong_efficiency_plot.jpeg")
  plt.clf()

  df = pd.read_csv("weak_scalability.csv")
  processes = df['#mpi processes'].values
  runtimes = df['total runtime'].values
  speedup = []
  for counter in range(1,6,2):
    speedup.append(runtimes[counter] / runtimes[counter -1])

  print(speedup)

  efficiency = [speedup[i]/processes[::2][i] for i in range(0,len(processes[::2]))]
  plt.scatter(processes[::2], efficiency, label='weak efficiency')
  plt.plot(processes[::2], efficiency)
  plt.plot(processes[::2], [1 for _ in processes[::2]], label='ideal weak efficiency')
  plt.ylabel('efficiency')
  plt.xlabel('#mpi processes')
  plt.legend(loc="lower left")
  plt.savefig("weak_efficiency_plot.jpeg")
  plt.clf()

  plt.scatter(processes[::2], speedup, label='weak scaling')
  plt.plot(processes[::2], speedup)
  plt.plot(processes[::2], processes[::2], label='linear speedup')
  plt.ylabel('speedup')
  plt.xlabel('#mpi processes')
  plt.legend(loc="upper left")
  plt.savefig("weak_scaling_plot.jpeg")
  plt.clf()
