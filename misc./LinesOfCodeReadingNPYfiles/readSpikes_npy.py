import numpy as np
import scipy.io
import pandas as pd


currentData = np.load('RatN_Day2_2019-10-11_03-58-54_spikes.npy', allow_pickle=True).item()

justInfo = currentData['info']

justTimes = dict()
justTimes['times'] = currentData['times']


pd.DataFrame(justInfo).to_csv("RatN_Day2_2019-10-11_03-58-54_spikesInfo.csv")
scipy.io.savemat('RatN_Day2_2019-10-11_03-58-54_spikesTimes.mat', justTimes)






