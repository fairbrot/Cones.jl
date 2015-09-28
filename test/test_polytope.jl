# Unit 2D square
A = [[1 0],
     [0 1],
     [-1 0],
     [0 -1]]
b = [0, 0, -1, -1]
points, rays, bidrays = find_generators(A,b)
P = Polytope(points, rays, bidrays)
mins, maxs = min_max_projections(P)
t_mins, t_maxs = Rational{Int}[0//1, 0//1], Rational{Int}[1//1, 1//1]
@test t_mins == mins
@test t_maxs == maxs

# Single point
A = [[1 1],
     [-1 -1],
     [1 -1],
     [-1 1]]
b = [1, -1, 0, 0]
points, rays, bidrays = find_generators(A,b)
P = Polytope(points, rays, bidrays)
mins, maxs = min_max_projections(P)
t_mins, t_maxs = Rational{Int}[1//2, 1//2], Rational{Int}[1//2, 1//2]
@test t_mins == mins
@test t_maxs == maxs


# The Fathomless Toblerone
A = [[1 0 0],
     [0 1 0],
     [-1 -1 0],
     [0 0 -1]]
b = [0, 0, -2, 0]
points, rays, bidrays = find_generators(A,b)
P = Polytope(points, rays, bidrays)
mins, maxs = min_max_projections(P, [eye(Rational{Int}, 3) [1, 1, 0]])
t_mins, t_maxs = Rational{Int}[0//1, 0//1, -1//0, 0//1], Rational{Int}[2//1, 2//1, 0//1, 2//1]
@test t_mins == mins
@test t_maxs == maxs

# Half-space {(x,y,z): z ≥ 1}
A = [[0 0 1]]
b = [1]
points, rays, bidrays = find_generators(A,b)
P = Polytope(points, rays, bidrays)
mins, maxs = min_max_projections(P)
t_mins, t_maxs = Rational{Int64}[-1//0,-1//0,1//1],Rational{Int64}[1//0,1//0,1//0]
@test t_mins == mins
@test t_maxs == maxs
