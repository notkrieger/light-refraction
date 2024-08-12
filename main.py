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
- y < 0 : not air, where we are interested in

WAVE:
- pick a point in space (y > 0) -- coming in from above
-  draw line to origin (centre of first boundary line)
- find perpendicular line, that is the plane wave moving toward the origin

MODEL:
- solve a system of simultaneous equations to figure refracted and transmitted light at each layer
- conserve momentum and energy between layers
"""

