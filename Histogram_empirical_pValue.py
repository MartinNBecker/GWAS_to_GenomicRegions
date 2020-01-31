# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 21:12:35 2020

@author: martbe5
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

permutation = pd.read_csv('TotalZscores_of_permutations.txt', header=None, names=['Zscore'])

compare = pd.read_csv('TotalZscore_of_tested_elements.txt', header=None, names=['Zscore'])


%matplotlib inline
plt.hist(permutation.Zscore , normed=True, bins=30)
plt.ylabel('Frequency');
plt.axvline(compare.Zscore[0], color = "r")

number_of_higher_perm_Z = len(permutation[permutation.Zscore >= compare.Zscore[0]])
number_of_perm = len(permutation.Zscore)
empir_p_value = number_of_higher_perm_Z / number_of_perm

