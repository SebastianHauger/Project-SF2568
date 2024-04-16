import matplotlib.pyplot as plt
import numpy as np
import transforms


def plot_data(data, compare=False, data2 = None):
    data = np.loadtxt(data, dtype=float)
    x = list(range(len(data)))
    plt.figure()
    plt.plot(x, data, label='computed')
    if compare:
        plt.plot(x, data2, label='correct')
        plt.plot(x, [x-y for x,y in zip(data, data2)], label='diff')
    plt.show()


if __name__=='__main__':
    compare = transforms.get_sol()
    plot_data('after_trans.txt', compare) 
    