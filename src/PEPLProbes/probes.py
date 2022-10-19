from PEPLProbes.utils.differentiate import derivative, SmoothDifferenceOperator
import numpy as np

def RPA(voltage, current, deriv_operator = SmoothDifferenceOperator()):
    '''
    Analyze an RPA trace

    Ported from code written by Leanne L. Su
    '''
    dI_dV = derivative(voltage, current, operator = deriv_operator)
    max_index = np.argmax(dI_dV)
    most_probable_voltage = voltage[max_index]

    return dI_dV, most_probable_voltage
