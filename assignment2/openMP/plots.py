import matplotlib.pyplot as plt
import pandas as pd

if __name__ == '__main__':
  df = pd.read_csv("strong_scalability.csv")
  threads = df['#threads'].values
  runtimes = df['runtime'].values
  speedup = [runtimes[0] / x for x in runtimes]

  plt.scatter(threads, speedup, label='strong scaling')
  plt.plot(threads, speedup)
  plt.ylabel('speedup')
  plt.xlabel('#threads')
  plt.legend(loc="upper left")
  plt.savefig("strong_scaling_plot.jpeg")
  plt.clf()

  df = pd.read_csv("weak_scalability.csv")
  threads = df['#threads'].values
  runtimes = df['runtime'].values
  speedup = [runtimes[0] / x for x in runtimes]

  plt.scatter(threads, speedup, label='weak scaling')
  plt.plot(threads, speedup)
  plt.ylabel('runtime [sec]')
  plt.xlabel('#threads')
  plt.legend(loc="upper left")
  plt.savefig("weak_scaling_plot.jpeg")
  plt.clf()
