from scipy.spatial import ConvexHull, convex_hull_plot_2d
from berakningsfysik.subsystem1.thermodynamics import entropy_at_endpoint
import numpy as np 
import matplotlib.pyplot as plt

print(entropy_at_endpoint(1))
x_values = [0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1]
G_mix_T1 = [1, 0.8, 0.9, 1.1, 1.2, 1.3, 1.21, 1.11, 2]
T_values = [0, 100, 200, 300]

points = np.array(list(zip(x_values, G_mix_T1)))

hull = ConvexHull(points)

x = np.arange(10)
print(x)
x.shape = (2,5)
print(x)
#plt.plot()
