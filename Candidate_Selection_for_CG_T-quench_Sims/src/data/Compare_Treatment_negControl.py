import pandas as pd
from scipy.stats import permutation_test, bootstrap
import numpy as np

################################################################################################################
# Load DataFrames from CSV files (replace 'file1.csv' and 'file2.csv' with your actual file paths)
df1 = pd.read_csv('data/Rand-False/EXP/Fract_High-OR_Cdists_system8_EXP_Rand-False_C.csv')
df2 = pd.read_csv('data/Rand-True/EXP/Fract_High-OR_Cdists_system8_EXP_Rand-True_C.csv')

# Select the column 'F' from each DataFrame
col_f1 = df1['F']
print(f'Treatment: {col_f1}')
col_f2 = df2['F']
print(f'NegControl: {col_f2}')

# Define the statistic function for the permutation test (difference in means)
def statistic(x, y):
    return x.mean() - y.mean()

# Perform the permutation test
result = permutation_test(
    (col_f1, col_f2),
    statistic,
    alternative='greater',
    n_resamples=10000,
    random_state=0
)

# Display the result
print(f'Treatment(mean) - NegControl(mean)')
print(f"Permutation test statistic: {result.statistic}")
print(f"P-value: {result.pvalue}")

################################################################################################################
#################### get the 95% confidence interval for the highest OR states
# Load DataFrames from CSV files (replace 'file1.csv' and 'file2.csv' with your actual file paths)
df1 = pd.read_csv('data/Rand-False/EXP/ranking_summary_system8_EXP_Rand-False_C.csv')
df1 = df1[df1['ranked_stateID'] == 0]
print(f'Treatment:\n{df1}')

df2 = pd.read_csv('data/Rand-True/EXP/ranking_summary_system8_EXP_Rand-True_C.csv')
df2 = df2[df2['ranked_stateID'] == 0]
print(f'Neg. Control:\n{df2}')

# Select the column 'F' from each DataFrame
col_f1 = df1['OR']
col_f2 = df2['OR']

# Define the statistic function for the permutation test (difference in means)
def statistic(x, y):
    return x.mean() - y.mean()

# Perform the permutation test
perm_result = permutation_test(
    (col_f1, col_f2),
    statistic,
    alternative='greater',
    n_resamples=10000,
    random_state=0
)

# Perform bootstrapping to get the 95% confidence intervals for each dataset
boot_result_f1 = bootstrap(
    (col_f1,),
    np.mean,
    confidence_level=0.95,
    n_resamples=10000,
    method='percentile',
    random_state=0
)

boot_result_f2 = bootstrap(
    (col_f2,),
    np.mean,
    confidence_level=0.95,
    n_resamples=10000,
    method='percentile',
    random_state=0
)

# Display results
print(f"Permutation test statistic: {perm_result.statistic}")
print(f"P-value: {perm_result.pvalue}")
print(f"95% Confidence interval for mean of Treatment: {np.mean(col_f1.values)} {boot_result_f1.confidence_interval}")
print(f"95% Confidence interval for mean of Neg. Control: {np.mean(col_f2.values)} {boot_result_f2.confidence_interval}")