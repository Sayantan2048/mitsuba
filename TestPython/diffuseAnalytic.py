import numpy as np
# Let's assume our shading point is at (0,0,0) and shading normal pointing at (0,0,1)
#2019-01-25 13:44:56 INFO  wrk0 [AnalyticDiffuseIntegrator] (-0.425173 0.905110 0.001966), (-0.733422 0.679766 0.003391), (-0.904923 0.425570 0.002123)
#2019-01-25 14:06:44 INFO  wrk3 [AnalyticDiffuseIntegrator] (-0.431031 0.902335 0.001963), (-0.739934 0.672671 0.003371), (-0.906886 0.421371 0.002111)

v1 = np.array([-0.431031, 0.902335, 0.001963])
v2 = np.array([-0.739934, 0.672671, 0.003371])
v3 = np.array([-0.906886, 0.421371, 0.002111])

#v1 = np.array([-0.425173, 0.905110, 0.001966])
#v2 = np.array([-0.733422, 0.679766, 0.003391])
#v3 = np.array([-0.904923, 0.425570, 0.002123])

#Normalize
v1 = v1 / np.sqrt(np.dot(v1,v1))
v2 = v2 / np.sqrt(np.dot(v2,v2))
v3 = v3 / np.sqrt(np.dot(v3,v3))

def integrateEdge(e1, e2):
    cosTheta = np.dot(e1, e2)
    theta = np.arccos(cosTheta)
    print("Theta:")
    print(theta)
    d = np.cross(e1, e2)[2]
    print("D:")
    print(d)
    t = theta / np.sin(theta)
    print("T:")
    print(t)
    return t * d

print(integrateEdge(v1, v2))
print(integrateEdge(v2, v3))
print(integrateEdge(v3, v1))
print(integrateEdge(v1, v2) + integrateEdge(v2, v3) + integrateEdge(v3, v1))