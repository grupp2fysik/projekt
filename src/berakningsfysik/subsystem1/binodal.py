from scipy.spatial import ConvexHull, convex_hull_plot_2d
#from berakningsfysik.subsystem1.thermodynamics import entropy_at_endpoint
import numpy as np 
import matplotlib.pyplot as plt

x_values = [0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1]
G_mix_T1 = [1, 0.8, 0.9, 1.1, 1.2, 1.3, 1.21, 1.11, 2]
T_values = [0, 100, 200, 300]

points = np.array(list(zip(x_values, G_mix_T1)))  # punkter i G-x-planet

hull = ConvexHull(points)  # det konvexa höljet till delta_G_mix(x)

plt.plot(points[:,0], points[:,1], 'o')
print(points)
print(hull.simplices)
print(hull.vertices)
# "simplices" är en vektor av index till par av
# punkter som utgör ändpunkter till höljets fasetter
for simplex in hull.simplices:
    plt.plot(points[simplex, 0], points[simplex, 1], 'k-') # hämta kolumn 0/1 för värdena som ingår i simplex

plt.title("Konvext hölje för delta_G_mix")
plt.xlabel("x")
plt.ylabel("delta_G_mix")
plt.savefig("convex_hull_G")


