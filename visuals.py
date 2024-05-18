import matplotlib.pyplot as plt
import numpy as np
import transforms


def plot_data():
    with open('timevals.txt', 'r') as f:
        lines = f.readlines()

    # Filter out the lines that start with 'n_processes'
    data_lines = [line for line in lines if not line.startswith('n_processes')]

    # Extract the floating point numbers
    data = [float(line.split(':')[1].strip()) for line in data_lines]
    parts = [[data[5*i + j] for i in range(int(len(data)/5))] for j in range(5)]
    total = [sum(parts[j][i] for j in range(len(parts))) for i in range(len(parts[0]))]
    # first_part = [data[5*i] for i in range(int(len(data)/5))]
    # second_part = [data[5*i + 1] for i in range(int(len(data)/5))]
    # third_part = [data[5*i + 2] for i in range(int(len(data)/5))]
    # fourth_part = [data[5*i + 3] for i in range(int(len(data)/5))]
    # fifth_part = [data[5*i + 4] for i in range(int(len(data)/5))]
    # total = [first_part[i] + second_part[i] + third_part[i] + fourth_part[i] + fifth_part[i] for i in range(len(first_part))]
    x = [i +1 for i in range(len(parts[0]))]
    print(sum(total))
    plt.figure()
    # Plot the data
    plt.plot(x, total)
    plt.show()
    plt.figure()
    labels = ["startup", "first shift", "DFT", "second shift", "idft"]
    for i,part in enumerate(parts):
        plt.plot(x, part, label=labels[i])
    plt.legend()
    plt.show()
    
    
def plot_bva():
    before = np.loadtxt("untouched.txt")
    after = np.loadtxt("transed_back.txt")
    plt.figure()
    plt.plot(before, linewidth=5.0)
    plt.plot(after)
    plt.show()


if __name__=='__main__':
    plot_bva()
    # plot_data()
    