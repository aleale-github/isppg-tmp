# Import Toolboxes and Input Data

import numpy as np


# Definition of the Payoff Function

def worst_off_autocallable_payoff(para_prod, risk_factor, barrier_factor,num_coupon_fix,num_sim):
    
    # Index Calculation
    rf1 = np.divide(risk_factor[:,:,1:],risk_factor[:,:,[0]])
    
    # Correction for Worst of
    rf = np.amin(rf1,axis=1)
    bf1 = barrier_factor / risk_factor[:,:,0]
    wo1 = np.min(bf1, axis=1)
    bf = wo1
        
    # Digital Effect
    pointer_digital_effect = para_prod[-(num_coupon_fix*3):-(num_coupon_fix*2), None].T        
    coupon_digital = para_prod[(num_coupon_fix*2):(num_coupon_fix*3)]
    digital_effect = (np.greater_equal(rf, para_prod[:num_coupon_fix]) * coupon_digital ) * pointer_digital_effect
    
    # Autocall Effect
    pointer_autocall_effect = para_prod[-(num_coupon_fix*2):-(num_coupon_fix), None].T
    
    # Pointer Matrix
    pointer_auto_call = np.greater_equal(rf, para_prod[num_coupon_fix:num_coupon_fix*2]).astype(int)
    pointer_auto_call[:, np.where(pointer_autocall_effect == 0)] = 0
    pointer_auto_call = np.cumsum(pointer_auto_call, axis=1)
    pointer_auto_call[np.where(pointer_auto_call > 1)] = 30000
    pointer_auto_call[np.where(pointer_auto_call[:,-1] == 0)] = 10000
    pointer_auto_call[np.where(pointer_auto_call == 10000)] = 3
    pointer_auto_call[np.where(pointer_auto_call == 20000)] = 2
    pointer_auto_call[np.where(pointer_auto_call == 30000)] = 0
    
    # Autocall Coupon T-
    auto_call_coupon = digital_effect
    auto_call_coupon[np.where(pointer_auto_call == 1)] += 1
    auto_call_coupon[np.where(pointer_auto_call == 0)] = 0
    
    # Autocall Coupon T
    pointer_auto_call_end = pointer_auto_call[:,-1]
    payoff_end = np.copy(bf)
    payoff_end[np.greater_equal(bf,para_prod[num_coupon_fix*2-1])] = para_prod[num_coupon_fix*3]
    payoff_end[np.greater_equal(bf,para_prod[num_coupon_fix-1])] = para_prod[num_coupon_fix*3+1]
    payoff_end[np.where(pointer_auto_call_end != 3)] = 0
    auto_call_coupon[:,-1] = payoff_end
    
    # Fix Coupon
    coupon_fix = np.repeat(para_prod[-(num_coupon_fix):][np.newaxis].T, num_sim, axis=1).T
    payoff = auto_call_coupon
    
    return payoff, coupon_fix
