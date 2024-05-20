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
    plt.title("Total time spent in a run of the better implementation")
    plt.ylabel("Time")
    plt.xlabel("#processes")
    plt.savefig("totTimeGood.pdf", dpi=200)
    plt.show()
    plt.figure()
    labels = ["Startup", "1. bitmirror", "FFT", "2. bitmirror", "IFFT"]
    for i,part in enumerate(parts):
        plt.plot(x, part, label=labels[i])
    plt.legend()
    plt.title("Decomposition of time spent in a run of the better implementation")
    plt.xlabel("#processes")
    plt.ylabel("time")
    plt.savefig("decTimeGood.pdf", dpi=200)
    plt.show()
    
    
def plot_bva():
    before = np.loadtxt("untouched_multiple_freq.txt")
    after = np.loadtxt("transed_back_multiple_freq.txt")
    plt.figure()
    plt.plot(before, linewidth=5.0, label="before transform")
    plt.plot(after, label="after both transforms")
    plt.title("Comparing initial data to data after forwards and backwards transform")
    plt.xlabel("Datapoint")
    plt.ylabel("Value")
    plt.legend()
    plt.savefig("multipleFreq.pdf", dpi=200)
    plt.show()

def get_efficiency():
    with open('timeValuesOldShift.txt', 'r') as f:
        lines = f.readlines()
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
    speedup = [total[0]/total[j] for j in range(len(total))]
    efficiency = [speedup[j]/(j+1) for j in range(len(speedup))]
    cost = [(j+1)*total[j] for j in range(len(total))]
    
    part = [sum(parts[j][i] for j in range(1,len(parts))) for i in range(len(parts[0]))]
    part_speedup = [part[0]/part[j] for j in range(len(part))]
    part_efficiency = [part_speedup[j]/(j+1) for j in range(len(part_speedup))]
    part_cost = [(j+1)*part[j] for j in range(len(total))] 
    
    # plt.figure()
    # plt.plot(x, speedup)
    # plt.title("Speedup of better code with startup")
    # plt.xlabel("#processes")
    # plt.ylabel("value")
    # plt.savefig("goodTotSpeedStartup.pdf", dpi = 200)
    # plt.show()
    # plt.figure()
    # plt.plot(x, efficiency)
    # plt.title("Efficiency of better code with startup")
    # plt.xlabel("#processes")
    # plt.ylabel("value")
    # plt.savefig("goodTotEFFStartup.pdf", dpi = 200)
    # plt.show()
    # plt.figure()
    # plt.plot(x, cost)
    # plt.title("Cost of better code with startup")
    # plt.xlabel("#processes")
    # plt.ylabel("value")
    # plt.savefig("goodTotCostStartup.pdf", dpi = 200)
    # plt.show()
    plt.figure()
    plt.plot(x, part_speedup)
    plt.title("Speedup of worse code without startup")
    plt.xlabel("#processes")
    plt.ylabel("value")
    plt.savefig("badTotSpeed.pdf", dpi = 200)
    plt.show()
    plt.figure()
    plt.plot(x, part_efficiency)
    plt.title("Efficiency of worse code without startup")
    plt.xlabel("#processes")
    plt.ylabel("value")
    plt.savefig("badTotEFF.pdf", dpi = 200)
    plt.show()
    plt.figure()
    plt.plot(x, part_cost)
    plt.title("Cost of worse code without startup")
    plt.xlabel("#processes")
    plt.ylabel("value")
    plt.savefig("badTotCost.pdf", dpi = 200)
    plt.show()
    

def plot_transformed():
    transformed = np.loadtxt("two_freq_after_shif.txt")
    plt.figure()
    plt.plot(transformed)
    plt.title("Trnasformation of the sum of two cosine waves")
    plt.xlabel("datapoint")
    plt.ylabel("value")
    plt.savefig("aTmultipleFreq.pdf", dpi=200)
    plt.show()
    


if __name__=='__main__':
    # plot_transformed()
    # plot_bva()
    # plot_data()
    get_efficiency()
    