import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import dct

def read_float(filename): 
    """read a textfile of floats and return a list of the values"""
    val = []
    with open(filename, 'r') as file:
        for line in file:
            val.append(float(line.strip()))
    return val

def get_test_data():
    seq = np.array([0.035844,
0.427558,
0.964480,
0.007856,
0.029386,
0.888409,
0.487235,
0.957974,
0.667670,
0.527987,
0.878195,
0.825866,
0.338036,
0.366950,
0.332365,
0.065925,

])
    dct_seq = dct(seq)
    print(dct_seq)

if __name__=='__main__':
    get_test_data()
    # untouched = read_float("untouched.txt")
    # doubletrans = read_float("transed_back.txt")
    # plt.figure()
    # plt.plot(untouched)
    # plt.show()
    # plt.figure()
    # plt.plot(doubletrans)
    # plt.show()
    