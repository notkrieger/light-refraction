import copy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

resolution = 1000
size = 10

boundary = resolution // 2
x = size * np.linspace(-1, 1, resolution)
z = size * np.linspace(-1, 1, resolution)


X0, Z0 = np.meshgrid(x, z)

# e's have to be bigger then 1, other like saying light travels faster in some other medium then vacuum
e0 = 1
m0 = 1
e1 = 1.2
sgn = np.sign(e1)
if sgn < 0:
    e1 *= -1 ### treat refractive index as positive for simplicty of solution
kappa = e1/e0


v0 = 1
v1 = e0*v0/e1

# Ei = ei*exp(-jk0z)
# Er = er*exp(jk0z)
# Et = et*exp(-jk1z)

# k0 = w(e0*m0)^0.5 = 2*pi/lam0
w = 4 ## controls how much
k0 = w * (e0*m0)**0.5
lam0 = 2*np.pi/k0
# k1 = w(e1*m0)^0.5 = kappa^0.5/k0 = 2*pi/lam1
k1 = kappa*k0 ### ks line up with not using kappa^1/2
lam1 = 2*np.pi/k1

initial_offset = [-3, -3]
if initial_offset[1] != 0:
    theta0 = np.arctan(initial_offset[0]/initial_offset[1])
else:
    theta0 = 2*np.pi
theta1 = np.arcsin(e0 / e1 * np.sin(theta0))


next_offset = [0, 0]
next_offset[0] = initial_offset[0] / (v0 * np.cos(theta0)) * v1 * np.cos(theta1)
next_offset[1] = initial_offset[1] / (v0 * np.sin(theta0)) * v1 * np.sin(theta1)

f0 = v0 / lam0

w0 = 2*np.pi*f0


# lengths
# distance w1 travels in layer 0 - width
#l1 = next_offset[0]
# true distance
l0 = 0
l1 = ((next_offset[0]-initial_offset[0])**2+(next_offset[1]-initial_offset[1])**2)**0.5

shift1 = k1*l0

# if error greater then critical angle or brewsters angle, when e1 > e0

# for any 2 boundaries need to solve
R1 = (1 - kappa**0.5)/(1 + kappa**0.5)
R2 = -R1
T12 = 1 + R1
T21 = 1 + R2

matrix = 1/T12*np.ones((2, 2), dtype='complex') # wave transmission matrix
matrix[0][1] *= R1
matrix[1][0] *= R1

phase_shift = np.zeros((2,2), dtype='complex')
phase_shift[0][0] = np.exp(1j*shift1)
phase_shift[1][1] = np.exp(-1j*shift1)


#apply phase shift backward

cb2 = np.asarray([1, 0], dtype='complex')

transfer_matrix = np.matmul(matrix, phase_shift)

cb = np.matmul(transfer_matrix, cb2)

Et0 = cb[0] * np.exp(-1j * (X0 * np.cos(theta0) + Z0 * np.sin(theta0)) * k0)
Er0 = cb[1] * np.exp(1j * (-X0 * np.cos(theta0) + Z0 * np.sin(theta0))*k0)
Et1 = cb2[0] * np.exp(-1j * (X0 * np.cos(theta1) + sgn * Z0 * np.sin(theta1)) * k1)
Er1 = cb2[1] * np.exp(1j * (X0 * np.cos(theta1) + Z0 * np.sin(theta1)) * k1) # last layer -- zero reflection

fig, ax = plt.subplots()
norm = plt.Normalize(-1, 1)
dt = 0.05


max0ind = 0
max0 = 0
max1ind = 0
max1 = 0
for i in range(resolution):
    if Et0[i, boundary] > max0:
        max0ind = i

for i in range(resolution):
    if Et1[i][boundary] > max1:
        max1ind = i

frame = ax.imshow(np.real(np.zeros_like(Et0)), norm=norm, extent=[-size, size, -size, size])

# shift whole array -- to make peaks line up at start of time
def shift(array1, n):
    arr = np.zeros_like(array1)
    for j in range(len(array1[0])):
        arr[j] = array1[(j-n) % resolution]
    return arr


Et0 = shift(Et0, max0ind)
Et1 = shift(Et1, max1ind)


def update(t):
    time0 = np.exp(1j*w0*t*dt)
    time1 = np.exp(1j * sgn * w0 * t * dt) # need separate time factor in case negative refraction
    E0 = np.real((Et0 + Er0)*time0)
    E1 = np.real((Et1 + Er1)*time1)

    # to make each field work in each layer
    remove0 = np.concatenate((np.ones_like(X0)[:, :boundary], np.zeros_like(X0[:, boundary:])), 1)
    remove1 = np.concatenate((np.zeros_like(X0)[:, :boundary], np.ones_like(X0[:, boundary:])), 1)

    superposition = np.real(E0 * remove0 + E1 * remove1)
    intensity = np.abs(superposition**2)
    #ax.imshow(np.real(Et0) * gausst0, norm=norm, extent=[-size, size, -size, size])

    #ax.pcolor(x, z, np.real(Et), norm=norm)
    frame.set_data(intensity)
    return [frame]


plt.title("incident light from " + str(e0) + " to " + str(sgn*e1))
# Create the animation object
ani = animation.FuncAnimation(fig, update, frames=70, blit=True)

# Save the animation to a file
ani.save('animation2.gif', fps=10)





