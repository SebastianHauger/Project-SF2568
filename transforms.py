import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import dct, idct, fft, ifft


def get_sol(transform='DCT'):
    transforms = ['DCT', 'DFT', 'IDCT', 'IDFT']
    if transform not in transforms:
        print("transform must be 'DCT', 'DFT', 'IDCT' or 'IDFT")
        return
    if transform==transforms[0]:
        data = np.loadtxt('untouched.txt', dtype=float)
        data = dct(data, type=3)
    elif transform==transforms[1]:
        data = np.loadtxt('untouched.txt', dtype=float)
        data = fft(data)
    elif transform==transforms[2]:
        data = np.loadtxt('after_transform.txt', dtype=float)
        data = idct(data, type=3)
    elif transform==transforms[3]:
        data = np.loadtxt('after_trans.txt', dtype=float)
        data = ifft(data)
    return data


if __name__ == '__main__':
    data = [0, 1, 2, 3, 4, 5, 6, 7]
    print(fft(data))
    
        