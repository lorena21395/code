import matplotlib.pyplot as plt
import numpy as np
plt.switch_backend('agg')

y = [-0.00673649,-0.0117174,-0.00195518,-0.00383003,0.00114767,0.00522412]
x = [0.001,0.01,0.1,1.,5,10]
err = [0.000471373,0.000555638,0.000514069,0.000203546,0.00291781,0.000577722]

plt.semilogx(x,y)
plt.errorbar(x,y,yerr = err)
plt.xlabel("log(Bg rms)")
plt.ylabel("Bias")
plt.show()
plt.savefig("noise_bias.png")
