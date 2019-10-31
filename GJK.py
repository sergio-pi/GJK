from numpy import cos, sin, pi, dot, cross, array, negative, array_equal, subtract, divide, transpose, linalg
from tkinter import Tk, ttk, Canvas

LARGE_FONT = ("Verdana", 12)
SMALL_FONT = ("Verdana", 9)
resolution = (800,600)
scale = 0.15

class Object():
    def __init__(self, points, **kwargs):
        self.localPoints = points
        self.worldPoints = list(points)
        self.translation = kwargs.get('translation', [0.0, 0.0])
        self.rotation = kwargs.get('rotation', [[1.0, 0.0], [0.0, 1.0]])
        self.angle = 0.0
        self.updateWorldSpace()
    def moveUp(self):
        self.translation[1] += 0.05
        self.updateWorldSpace()
    def moveLeft(self):
        self.translation[0] -= 0.05
        self.updateWorldSpace()
    def moveDown(self):
        self.translation[1] -= 0.05
        self.updateWorldSpace()
    def moveRight(self):
        self.translation[0] += 0.05
        self.updateWorldSpace()
    def rotateObject(self, theta):
        self.angle = (self.angle + theta) % (2*pi)
        matrix = rotationMatrix(self.angle)
        self.rotation = matrix
        self.updateWorldSpace()
    def updateWorldSpace(self):
        self.worldPoints = [vecMatrix(p, self.rotation) for p in  self.localPoints]         # Rotate
        self.worldPoints = [array(p) + array(self.translation) for p in self.worldPoints]   # Translate
    def localToWorld(self, v):
        v = vecMatrix(v, self.rotation)         # Rotate
        v = array(v) + array(self.translation)  # Translate
        return v
    def worldToLocal(self, v):
        v = vecMatrix(v, transpose(self.rotation))  # Inverse Rotation
        v = array(v) - array(self.translation)      # Inverse Translation
        return v
    def localToWorldVec(self, v):
        return vecMatrix(v, self.rotation)             # Rotate
    def worldToLocalVec(self, v):
        return vecMatrix(v, transpose(self.rotation))  # Inverse Rotation

class RegularPolygon(Object):
    def __init__(self, n, **kwargs):
        if n == 0: points = []
        angle = 2*pi/n
        points = [[cos(angle * t), sin(angle * t)] for t in range(n)]
        Object.__init__(self, points, **kwargs)

class Polygon(Object):
    def __init__(self, points, **kwargs):
        Object.__init__(self, points, **kwargs)

class MinkowskiDifference(Object):
    def __init__(self, A, B, **kwargs):
        points = []
        for i in range(30):
            newp = csoSupport(A, B, [cos(2*pi/30 * i), sin(2*pi/30 * i)])
            duplicate = False
            for s in points:
                if (s == newp).all(): duplicate = True
            if not duplicate: points.append(newp)
            
        Object.__init__(self, points, **kwargs)

class SceneObjects():
    def __init__(self):
        self.objects = []
    def addObject(self, *o):
        for i in o:
            self.objects.append(i)

def toCanvas(points):
    aspectRatio = float(resolution[1]) / float(resolution[0])
    return [[resolution[0] * (0.5 + x * scale * aspectRatio),
            resolution[1] - (resolution[1] * (0.5 + y * scale))] for (x, y) in points]

def vecMatrix(v, m):
    transformed = []
    for vdim in range(len(v)):
        dotp = sum([a*b for a, b in zip(v, [i[vdim] for i in m])])
        transformed.append(dotp)
    return transformed

def normalize(v):
    norm = linalg.norm(array(v))
    nv = [0.0, 0.0]
    if norm != 0:
        nv = divide(v, norm)
    return nv

def rotationMatrix(theta):
    return ((cos(theta), sin(theta)), (-sin(theta), cos(theta)))

def support(v, s):
    maxdot = 0.0
    maxpoint = None
    for i in s:
        currentDot = dot(v, i)
        if currentDot > maxdot:
            maxdot = currentDot
            maxpoint = i
    return maxpoint

def barycentricEdge(A, B, Q):
    # Compute barycentric coordinates u and v
    n = subtract(B,A)
    v = divide(dot(subtract(Q,A), n), dot(n, n))
    u = 1.0 - v
    return u,v

def triangleArea(A, B, C):
    return 0.5 * cross(subtract(B,A), subtract(C,A))

def closestPointOnEdgeToPoint(A, B, Q):
    u,v = barycentricEdge(A, B, Q)
    
    if u <= 0:
        return B
    elif v <= 0:
        return A
    else:
        return u*A + v*B

def closestPointOnTriangleToPoint(A, B, C, Q):
    # Compute barycentric coordinates of all the edges
    uAB, vAB = barycentricEdge(A, B, Q)
    uBC, vBC = barycentricEdge(B, C, Q)
    uCA, vCA = barycentricEdge(C, A, Q)
    # Test the vertex regions
    if vAB <= 0 and uCA <= 0: return A, 0  # Region A
    if uAB <= 0 and vBC <= 0: return B, 1  # Region B
    if uBC <= 0 and vCA <= 0: return C, 2  # Region C
    # Compute barycentric coordinates of the triangle
    areaABC = triangleArea(A,B,C)
    #if areaABC == 0.0: print "Div by 0" TODO: handle triangles with zero area that trigger infinite loop
    uABC = triangleArea(Q,B,C) / areaABC
    vABC = triangleArea(Q,C,A) / areaABC
    wABC = triangleArea(Q,A,B) / areaABC    # can be replaced with 1 - uABC - wABC
    # Test the edge regions
    if uBC > 0 and vBC > 0 and uABC <= 0: return uBC * B + vBC * C, 3  # Region BC
    if uCA > 0 and vCA > 0 and vABC <= 0: return uCA * C + vCA * A, 4  # Region CA
    if uAB > 0 and vAB > 0 and wABC <= 0: return uAB * A + vAB * B, 5  # Region AB
    # Region ABC (uABC > 0 and vABC > 0 and wABC > 0)
    return Q, 6

def perpToEdgePointingTowards(A, B, Q):
    d = cross(subtract(A, B), [0.0, 0.0, 1.0])
    d = d[:-1]
    if dot(subtract(Q, A), d) < 0:
        d = negative(d)
    return d

def GJK(self, A, B):
    simplex = []
    d = [0.7071, 0.7071]    # Initial arbitrary direction
    P = []                  # Point of minimum norm in simplex
    origin = [0.0, 0.0]
    r = -1                  # Voronoi region of the triangle that tells how to mutate the simplex
    simplex.append(csoSupport(A, B, d))
    while(1):
        simplexLen = len(simplex)
        if simplexLen == 1:
            P = simplex[0]
            d = negative(P)
        elif simplexLen == 2:
            P = closestPointOnEdgeToPoint(simplex[0], simplex[1], origin)
            d = perpToEdgePointingTowards(simplex[0], simplex[1], origin)
        elif simplexLen == 3:
            P, r = closestPointOnTriangleToPoint(simplex[0], simplex[1], simplex[2], origin)
        visualization(self, P, d, simplex)
        if array_equal(P, origin): return 0.0   # Exit if the origin is inside the triangle
        # Mutate simplex and update the next search d
        if simplexLen == 3:
            if r < 3:   # The voronoi region is a vertex, use such vertex as the new simplex
                simplex = [simplex[r]]
                d = negative(simplex[0])
            else:       # The voronoi region is an edge, get rid of one vertex from the simplex
                del simplex[r - 3]
                d = perpToEdgePointingTowards(simplex[0], simplex[1], origin)
        V = csoSupport(A, B, d)
        if any((V == x).all() for x in simplex): return linalg.norm(array(P))
        simplex.append(V)

def csoSupport(colliderA, colliderB, dir):
    # Convert search dir to world space
    localDirA = colliderA.worldToLocalVec(dir)
    localDirB = colliderB.worldToLocalVec(negative(dir))
    # Compute support points in local space
    supportA = support(localDirA, colliderA.localPoints)
    supportB = support(localDirB, colliderB.localPoints)
    # Convert support points to world space
    supportA = colliderA.localToWorld(supportA)
    supportB = colliderB.localToWorld(supportB)
    # Compute CSO support point
    return subtract(supportA, supportB)

def visualization(self, P, d, simplex):
    self.searchDirection.worldPoints = [ P, P + d]
    self.simplex.worldPoints = simplex

class MainApp(Tk):
    def __init__(self, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)
        self.wm_title("Boolean GJK Algorithm")
        self.minsize(width=resolution[0], height=resolution[1])
        self.tcanvas = Canvas(self, background = "black")
        self.tcanvas.pack(fill="both", expand=True)
        self.bind('<q>', self.rotateInput)
        self.bind('<e>', self.rotateInput)
        self.bind('<w>', self.wasdInput)
        self.bind('<a>', self.wasdInput)
        self.bind('<s>', self.wasdInput)
        self.bind('<d>', self.wasdInput)

        originLine1 = Object([[-0.25, 0.0], [0.25, 0.0]])
        originLine2 = Object([[0.0, -0.25], [0.0, 0.25]])

        self.collision = False
        self.triangle = RegularPolygon(3)
        self.box = Polygon([[1.0, -1.0],[1.0, 1.0],[-1.0, 1.0],[-1.0, -1.0]], translation = [3.0, 0.0])
        self.myscene = SceneObjects()
        self.myscene.addObject(self.triangle, self.box, originLine1, originLine2)

        self.searchDirection = Polygon([[0.0, 0.0], [0.0, 0.0]])
        self.simplex = Polygon([[0.0, 0.0], [0.0, 0.0]])
        self.minkowski = MinkowskiDifference(self.triangle, self.box)

        self.myscene.addObject(self.minkowski, self.searchDirection, self.simplex)

        self.drawScene(self.myscene)

    def rotateInput(self, event):
        theta = 0.05
        if event.keycode == 81:
            self.triangle.rotateObject(theta)
        elif event.keycode == 69:
            self.triangle.rotateObject(-theta)
        self.checkCollision(self.triangle, self.box)
        self.drawScene(self.myscene)

    def wasdInput(self, event):
        if event.keycode == 87:
            self.triangle.moveUp()
        if event.keycode == 65:
            self.triangle.moveLeft()
        if event.keycode == 83:
            self.triangle.moveDown()
        if event.keycode == 68:
            self.triangle.moveRight()
        self.checkCollision(self.triangle, self.box)
        self.drawScene(self.myscene)

    def checkCollision(self, A , B):
        self.collision = True if GJK(self, A, B) <= 0.0 else False

    def drawScene(self, scene):
        self.minkowski.__init__(self.triangle, self.box)
        self.tcanvas.delete('all')
        for obj in scene.objects:
            canvasPoints = toCanvas(obj.worldPoints)   # ready to print geometry in the canvas space
            for i,v in enumerate(canvasPoints):
                if i == len(canvasPoints) - 1: i = 0
                else: i += 1
                self.tcanvas.create_line("{0} {1} {2} {3}".format(
                v[0], v[1], canvasPoints[i][0], canvasPoints[i][1]
                    ), fill= "#FF1111" if self.collision else "#ACACAC")
app = MainApp()
app.mainloop()