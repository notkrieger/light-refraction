# refraction simulation
"""
GOALS:
- make simple 2d refraction for single source plane wave going between two materials
- after this there are multiple options:
     - make for negative refractive index materials
     - make for multiple layers of materials
     - make for multiple waves from different points
- ideally we want multiple of these options all cooperating together

SETUP:
- y > 0 : in the air
- y = 0 : first boundary
- y < 0 : other material/s, where we are interested in

WAVE:
- pick a point in space (y > 0) -- coming in from above
-  draw line to origin (centre of first boundary line)
- find perpendicular line, that is the plane wave moving toward the origin

MODEL:
- solve a system of simultaneous equations to figure refracted and transmitted light at each layer
- conserve momentum and energy between layers
"""


### https://www.dev-mind.blog/simulating-light/
### https://github.com/ngmsoftware/FD-TD

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


fig, ax = plt.subplots()

size = 5
x = size*np.linspace(-1, 1, 200)
y = size*np.linspace(-1, 1, 200)

X, Y = np.meshgrid(x, y)


k = [-1, -1]
#amp = np.exp(-np.linspace(-np.pi, np.pi, 100)**2)

F = np.exp(1j*np.pi*(k[0]*X+k[1]*Y))

E = 0
H = 0

dt = 0.1
for t in range(10):
    data = (np.real(F * np.exp(1j*np.pi*2*dt*t)))

    ax.imshow(data)
    #imagesc(x, z, reshape(sum(real(F * exp(1i * 2 * pi * dt * ii)), 2), size(X)));
    plt.pause(0.000001)
plt.show()
