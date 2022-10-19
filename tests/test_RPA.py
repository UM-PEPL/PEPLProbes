from PEPLProbes.probes import RPA
from PEPLProbes.utils.differentiate import SmoothDifferenceOperator, CentralDifferenceOperator
import numpy as np
from scipy.special import erf
import math

def rpa_helper(deriv_operator):
    # Sweep voltage from 0 to 300 V
    voltage = np.linspace(0, 300, 301)

    # dI_dV has a mean of mu and a standard deviation of sigma
    sigma = 100.0
    mu = 150.0

    # Define sample current trace as a sigmoid (error function)
    current = 0.5 * (1 + erf((voltage - mu)/sigma/np.sqrt(2)))

    # Derivative of erf is a Gaussian
    expected_dI_dV = np.exp(-0.5 * (voltage - mu)**2 / sigma**2) / sigma / np.sqrt(2 * math.pi)

    # Compute numerical derivative of trace
    dI_dV, most_probable_voltage = RPA(voltage, current, deriv_operator)

    # Check that most probable voltage is correct (mean of Gaussian)
    assert most_probable_voltage == mu, "Most probable voltage is incorrect"

    # Check that derivative is correct
    assert np.allclose(dI_dV, expected_dI_dV, atol = min(expected_dI_dV)), "Derivative of trace is incorrect"

def test_rpa_smooth():
    rpa_helper(SmoothDifferenceOperator())

def test_rpa_central():
    rpa_helper(CentralDifferenceOperator())
