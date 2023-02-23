import numpy as np
import scipy.io
import pandas as pd


currentData = np.load('RatN_Day2_2019-10-11_03-58-54.tracks.laps.npy', allow_pickle=True).item()

scipy.io.savemat('RatN_Day2_2019-10-11_03-58-54.tracks.laps.mat', currentData)






