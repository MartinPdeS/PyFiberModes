"""
Coupled-mode power exchange
===========================

Propagate two lossless modes coupled by a Hermitian matrix.
"""

import matplotlib.pyplot as plt
import numpy as np
from PyFiberModes.units import ureg

from PyFiberModes.propagation import CoupledModeSystem

z = np.linspace(0, 40, 400) * ureg.meter
system = CoupledModeSystem(betas=np.array([0, 0]) / ureg.meter, coupling=np.array([[0, 0.12], [0.12, 0]]) / ureg.meter)
result = system.propagate(initial=[1, 0], z=z)

figure, axis = plt.subplots()
axis.plot(z.to("meter").magnitude, result.powers[:, 0], label="mode 1")
axis.plot(z.to("meter").magnitude, result.powers[:, 1], label="mode 2")
axis.set(xlabel="Propagation distance", ylabel="Normalized power", ylim=(0, 1))
axis.legend()
plt.show()
