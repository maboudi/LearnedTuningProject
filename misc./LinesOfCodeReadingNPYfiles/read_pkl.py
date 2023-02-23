import numpy as np
import scipy.io
import pandas as pd


currentData = pd.read_pickle('RatN_Day2_2019-10-11_03-58-54.states.pkl')

pd.DataFrame(currentData).to_csv('RatN_Day2_2019-10-11_03-58-54.states.csv')


# scipy.io.savemat('RatN_Day2_2019-10-11_03-58-54.states.mat', currentData)






