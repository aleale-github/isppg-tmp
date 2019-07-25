# Import Toolboxes and Input Data
import sys
sys.path.append("../")
import numpy as np
import pandas as pd
import worst_off_autocallable as woa
import scipy.io as spio

def test_worst_off_autocallable():

    mat = spio.loadmat('../data/WorstOffAutocallablePayoff_inputs.mat', squeeze_me=True)

    # Function Execution
    Payoff, CouponFix = woa.worst_off_autocallable_payoff(mat['ParaProd'], mat["RiskFactor"], mat["BarrierFactor"], mat["NumCouponFix"], mat["NumSim"])

    # Load Matlab Results
    MatPayoff = pd.read_csv('../data/Payoff.csv', header=None).to_numpy()
    MatCouponFix = pd.read_csv('../data/CouponFixed.csv', header=None).to_numpy()

    assert np.isclose(Payoff,MatPayoff).all()
    assert np.isclose(CouponFix, MatCouponFix).all()

    # Payoff Results
    # results = pd.DataFrame({"Payoff Python": Payoff[:,-1], "Payoff Matlab": MatPayoff[:,-1]})
    # print(results)

