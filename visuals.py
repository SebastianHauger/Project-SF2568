import matplotlib.pyplot as plt


def read_float(filename): 
    """read a textfile of floats and return a list of the values"""
    val = []
    with open(filename, 'r') as file:
        for line in file:
            val.append(float(line.strip()))
    return val

if __name__=='__main__':
    untouched = read_float("untouched.txt")
    doubletrans = read_float("transed_back.txt")
    plt.figure()
    plt.plot(untouched)
    plt.show()
    plt.figure()
    plt.plot(doubletrans)
    plt.show()
    