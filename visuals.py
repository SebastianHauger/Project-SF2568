import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import dct, idct

def read_float(filename): 
    """read a textfile of floats and return a list of the values"""
    val = []
    with open(filename, 'r') as file:
        for line in file:
            val.append(float(line.strip()))
    return val

def get_test_data():
    seq = np.array([0.547858,
0.848343,
0.093534,
0.023399,
0.264044,
0.780526,
0.303612,
0.808418,
0.076458,
0.037411,
0.773003,
0.859976,
0.623247,
0.908122,
0.806442,
0.877579

])
    dct_seq = dct(seq)
    seq = idct(dct_seq)
    print("dct seq", dct_seq)
    print("inverse dct seq", seq)

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
    