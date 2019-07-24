# Import Toolboxes and Input Data

import scipy.io as spio
import numpy as np
import pandas as pd
import worst_off_autocallable as woa

mat = spio.loadmat('WorstOffAutocallablePayoff_inputs.mat', squeeze_me=True)

# Function Execution

Payoff, CouponFix = woa.worst_off_autocallable_payoff(mat['ParaProd'], mat["RiskFactor"], mat["BarrierFactor"], mat["NumCouponFix"], mat["NumSim"])
Payoff = pd.DataFrame(Payoff)
CouponFix = pd.DataFrame(CouponFix)

# Load Matlab Results

MatPayoff = pd.read_csv('Payoff.csv', header=None)
MatCouponFix = pd.read_csv('CouponFixed.csv', header=None)

# Differences between Python and Matlab

payoff_diff = ((Payoff - MatPayoff) / MatPayoff).round(15)
payoff_diff.fillna(0, inplace=True)
coupon_fix_diff = ((CouponFix - MatCouponFix) / MatCouponFix).round(15)
coupon_fix_diff.fillna(0, inplace=True)

# Payoff Results

print('\n' + 'Payoff Results' + '\n')
results = pd.DataFrame(np.transpose([Payoff.iloc[:,-1], MatPayoff.iloc[:,-1], payoff_diff.iloc[:,-1].abs()]), columns=["Payoff Python", "Payoff Matlab", "Differences"])
print(results)

# Payoff differences Min and Max

print('\n' + 'Payoff differences Min and Max' + '\n')
payoff_min_max = pd.DataFrame([payoff_diff.min(), payoff_diff.max()])
print(payoff_min_max)

# CouponFix differences Min and Max

print('\n' + 'CouponFix differences Min and Max' + '\n')
coupon_fix_min_max = pd.DataFrame([coupon_fix_diff.min(), coupon_fix_diff.max()])
print(coupon_fix_min_max)
