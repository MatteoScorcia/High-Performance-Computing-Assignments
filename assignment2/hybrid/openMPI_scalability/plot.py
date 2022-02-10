import matplotlib.pyplot as plt
import pandas as pd

if __name__ == '__main__':
  df = pd.read_csv("strong_scalability.csv")
  threads = df['#threads'].values
  runtimes = df['total runtime'].values
  speedup = [runtimes[0] / x for x in runtimes]

  print(speedup)

  plt.scatter(threads, speedup, label='strong scaling')
  plt.plot(threads, speedup)
  plt.plot(threads, threads, label='linear speedup')
  plt.ylabel('speedup')
  plt.xlabel('#threads')
  plt.legend(loc="upper left")
  plt.savefig("strong_scaling_plot.jpeg")
  plt.clf()

  efficiency = [speedup[i]/threads[i] for i in range(0,len(threads))]
  plt.scatter(threads, efficiency, label='strong efficiency')
  plt.plot(threads, efficiency)
  plt.plot(threads, [1 for _ in threads], label='ideal strong efficiency')
  plt.ylabel('efficiency')
  plt.xlabel('#threads')
  plt.legend(loc="lower left")
  plt.savefig("strong_efficiency_plot.jpeg")
  plt.clf()

  df = pd.read_csv("weak_scalability.csv")
  threads = df['#threads'].values
  runtimes = df['total runtime'].values
  speedup = []
  for counter in range(1,10,2):
    speedup.append(runtimes[counter] / runtimes[counter -1])

  print(speedup)

  efficiency = [speedup[i]/threads[::2][i] for i in range(0,len(threads[::2]))]
  plt.scatter(threads[::2], efficiency, label='weak efficiency')
  plt.plot(threads[::2], efficiency)
  plt.plot(threads[::2], [1 for _ in threads[::2]], label='ideal weak efficiency')
  plt.ylabel('efficiency')
  plt.xlabel('#threads')
  plt.legend(loc="lower left")
  plt.savefig("weak_efficiency_plot.jpeg")
  plt.clf()

  plt.scatter(threads[::2], speedup, label='weak scaling')
  plt.plot(threads[::2], speedup)
  plt.plot(threads[::2], threads[::2], label='linear speedup')
  plt.ylabel('speedup')
  plt.xlabel('#threads')
  plt.legend(loc="upper left")
  plt.savefig("weak_scaling_plot.jpeg")
  plt.clf()
