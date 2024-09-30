import copy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

resolution = 1000 # 1001?, worry about this later
size = 10

x = size * np.linspace(-1, 1, resolution)
z = size * np.linspace(-1, 1, resolution)

X, Z = np.meshgrid(x, z)

es = [1, 1.2, 3]
lens = [2, 2, 1] # lens of layers, to be normalised to size
if len(lens) != len(es):
    print('layer lengths and layer permeativities are different lengths')

num_layers = len(es)

# for ease make sum of lens a factor of resolution
# 1, 2, 4, 5, 8, 10, ...
lens /= np.sqrt(np.sum(np.asarray(lens))**2) ## normalise array

sgns = np.sign(es) # represents a change in polarity across a boundary

es = np.abs(es) # make es to remove problems with negatives
ms = np.ones_like(es) # permeabilities

boundaries = np.zeros(len(lens)-1) # defined with respect to resolution, not size
kappas = np.ones_like(boundaries) # boundary constants??

thetas = np.zeros_like(es)
thetas[0] = np.pi / 4 # initial angle of incidence

## solve thetas
for ii in range(1, num_layers):
    thetas[ii] = np.arcsin(es[ii-1] / es[ii] * np.sin(thetas[ii-1]))

# calculate boundary locations
for ii in range(len(lens)-1):
    kappas[ii] = es[ii+1]/es[ii]
    if ii == 0:
        boundaries[ii] = lens[ii] * resolution
    else:
        boundaries[ii] = boundaries[ii-1] + lens[ii] * resolution

w = 3.5 # angular frequency
ks = np.ones_like(es) # wavenumber for each layer
# calculate wave number in each layer
for ii in range(len(es)):
    if ii == 0:
        ks[ii] = w * (es[ii]*ms[ii])**0.5
    else:
        ks[ii] = ks[ii-1] * kappas[ii-1]



# check for bugs when adding more layers
ls = np.zeros(len(boundaries)) # lengths of layers to calculate phase shift, this is in size, not res
for ii in range(1, ls.size):
    ls[ii] = x[int(boundaries[ii])] - x[int(boundaries[ii-1])]

ts = np.ones(num_layers, dtype='complex') # transmitted at each layer
rs = np.zeros(num_layers, dtype='complex') # reflected at each layer


Et = np.zeros((es.size, resolution, resolution), dtype='complex') # construct transmission field for each layer
Er = np.zeros_like(Et, dtype='complex')

previous_transfer_matrix = np.ones((2,2), dtype='complex')
print(ls)
print(boundaries)
for ii in range(boundaries.size):

    R1 = (1 - kappas[num_layers - ii - 2] ** 0.5) / (1 + kappas[num_layers - ii - 2] ** 0.5)
    T12 = 1 + R1

    rt_matrix = 1/T12*np.ones((2, 2), dtype='complex')
    rt_matrix[0][1] *= R1
    rt_matrix[1][0] *= R1

    phase_shift = np.zeros((2, 2), dtype='complex')
    phase_shift[0][0] = np.exp(1j * ls[ii])
    phase_shift[1][1] = np.exp(-1j * ls[ii])

    #transfer_matrix = previous_transfer_matrix @ rt_matrix @ phase_shift
    transfer_matrix = np.matmul(rt_matrix, phase_shift)

    rt = np.asarray([ts[num_layers - ii - 1], rs[num_layers - ii - 1]])
    rt2 = np.matmul(transfer_matrix, rt)
    ts[num_layers - ii - 2] = rt2[0]
    rs[num_layers - ii - 2] = rt2[1]

    if ii != boundaries.size - 1:
        previous_transfer_matrix = copy.deepcopy(transfer_matrix)

print(rs)
print(ts)
t0 = np.abs(ts[0])
ts /= t0 #normalise with respect to inital wave
rs /= t0


for ii in range(num_layers):
    # adjust sign so its relative to previous boundary,
    # ie. has there been a switch in sign of refractive index accross the boundary
    Et[ii] = ts[ii] * np.exp(-1j * (X * np.cos(thetas[ii]) + sgns[ii] * Z * np.sin(thetas[ii])) * ks[ii])
    Er[ii] = rs[ii] * np.exp(1j * (-X * np.cos(thetas[ii]) + sgns[ii] * Z * np.sin(thetas[ii])) * ks[ii])

remove = np.zeros_like(Et)
prev_boundary = 0
for ii in range(num_layers):
    # solve so remove only makes ones in boundary for that layer
    if ii == num_layers-1:
        next_boundary = resolution
    else:
        next_boundary = int(boundaries[ii])

    for jj in range(prev_boundary, next_boundary):
        remove[ii, :, jj] += 1

    if ii != num_layers-1:
        prev_boundary = int(boundaries[ii])



# make animation
fig, ax = plt.subplots()
norm = plt.Normalize(-1, 1)
dt = 0.05

frame = ax.imshow(np.real(np.zeros_like(Et[0])), norm=norm, extent=[-size, size, -size, size])

def update(t):
    time = np.exp(1j*w*sgns*t*dt)

    E = np.zeros_like(Et[0])
    for ii in range(num_layers):
        E = E + (Et[ii] + Er[ii]) * time[ii] * remove[ii]

    superposition = np.real(E)
    intensity = np.abs(superposition**2)

    frame.set_data(superposition)
    return [frame]


# Create the animation object
ani = animation.FuncAnimation(fig, update, frames=100, blit=True)

# Save the animation to a file
ani.save('multiLayer.gif', fps=10)

