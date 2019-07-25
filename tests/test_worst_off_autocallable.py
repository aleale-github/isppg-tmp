# Import Toolboxes and Input Data
import sys
sys.path.append("../")
import numpy as np
import pandas as pd
import isppg.payoffs.worst_off_autocallable as woa
import scipy.io as spio
import os
datadir = os.path.abspath(os.path.dirname(woa.__file__)+"/../../data/")

def test_worst_off_autocallable():
    print(datadir)

    mat = spio.loadmat('{}/WorstOffAutocallablePayoff_inputs.mat'.format(datadir), squeeze_me=True)

    # Function Execution
    Payoff, CouponFix = woa.worst_off_autocallable_payoff(mat['ParaProd'], mat["RiskFactor"], mat["BarrierFactor"], mat["NumCouponFix"], mat["NumSim"])

    # Load Matlab Results
    MatPayoff = pd.read_csv('{}/Payoff.csv'.format(datadir), header=None).to_numpy()
    MatCouponFix = pd.read_csv('{}/CouponFixed.csv'.format(datadir), header=None).to_numpy()

    assert np.isclose(Payoff,MatPayoff).all()
    assert np.isclose(CouponFix, MatCouponFix).all()

    # Payoff Results
    results = pd.DataFrame({"Payoff Python": Payoff[:,-1], "Payoff Matlab": MatPayoff[:,-1]})
    print(results)

if __name__ == "__main__":
    test_worst_off_autocallable()
