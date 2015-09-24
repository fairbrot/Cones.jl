# Unit 2D square
A = [[1 0],
     [0 1],
     [-1 0],
     [0 -1]]
b = [0, 0, -1, -1]
points, rays, bidrays = find_generators(A,b)

# Single point
A = [[1 1],
     [-1 -1],
     [1 -1],
     [-1 1]]
b = [1, -1, 0, 0]
points, rays, bidrays = find_generators(A,b)

# Toblerone
A = [[1 0 0],
     [0 1 0],
     [-1 -1 0]]
b = [0, 0, -2]
points, rays, bidrays = find_generators(A,b)
