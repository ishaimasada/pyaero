import numpy


class Material:
    def __init__(self, youngs_modulus, shear_modulus, density, poisson_ratio):
        self.E = youngs_modulus
        self.G = shear_modulus
        self.rho = density
        self.nu = poisson_ratio

class ShaftElement:
    """
    Two-node rotor shaft element (Timoshenko / Euler-Bernoulli) along the z axis.

    DOFs per node: [u, v, theta_x, theta_y]  (u: x-disp, v: y-disp, theta_x/theta_y: rotations about x and y, right-hand rule).
    Slopes:  theta_y = +du/dz  (x-z plane),  theta_x = -dv/dz  (y-z plane).
    Element DOF order: [u1, v1, tx1, ty1, u2, v2, tx2, ty2]  (8 x 8 matrices).

    Flags:
        shear_effects  : True -> Timoshenko (Phi from kappa*G*A), False -> Euler-Bernoulli (Phi = 0)
        rotary_inertia : True -> add rotary inertia matrix to the mass matrix
        gyroscopic     : True -> gyroscopic_matrix() is non-zero (multiply by spin speed Omega)
    """

    # indices of the 4 DOFs of each bending plane in the 8-DOF element vector
    _IDX_A = [0, 3, 4, 7]   # x-z plane: u1, ty1, u2, ty2
    _IDX_B = [1, 2, 5, 6]   # y-z plane: v1, tx1, v2, tx2
    # sign flip: slope in y-z plane is -theta_x
    _T = numpy.diag([1.0, -1.0, 1.0, -1.0])

    def __init__(self, length, inner_diameter, outer_diameter, material:Material, rotary_inertia=False, shear_effects=False, gyroscopic=False):
        self.L = length
        self.ID = inner_diameter
        self.OD = outer_diameter
        self.material = material
        self.rotary_inertia = rotary_inertia
        self.shear_effects = shear_effects
        self.gyroscopic = gyroscopic

    # Section properties
    @property
    def area(self):
        """Cross-section area  A = pi/4 (OD^2 - ID^2)"""
        return numpy.pi / 4.0 * (self.OD**2 - self.ID**2)

    @property
    def inertia(self):
        """Second moment of area  I = pi/64 (OD^4 - ID^4)"""
        return numpy.pi / 64.0 * (self.OD**4 - self.ID**4)

    @property
    def polar_inertia(self):
        """Polar second moment of area  J = 2 I"""
        return 2.0 * self.inertia

    @property
    def shear_coefficient(self):
        """Cowper's shear correction factor for a hollow circular section"""
        nu = self.material.nu
        m = self.ID / self.OD
        a = (1.0 + m**2)**2
        return 6.0 * (1.0 + nu) * a / ((7.0 + 6.0 * nu) * a + (20.0 + 12.0 * nu) * m**2)

    @property
    def phi(self):
        """Shear parameter  Phi = 12 E I / (kappa G A L^2)  (0 for Euler-Bernoulli)"""
        if not self.shear_effects:
            return 0.0
        return 12.0 * self.material.E * self.inertia / ( self.shear_coefficient * self.material.G * self.area * self.L**2)

    # Single-plane (4x4) matrices, DOF order [w1, slope1, w2, slope2]
    def _stiffness_plane(self):
        E, I, L, p = self.material.E, self.inertia, self.L, self.phi
        c = E * I / (L**3 * (1.0 + p))
        return c * numpy.array([
            [ 12.0,      6.0*L,         -12.0,      6.0*L],
            [ 6.0*L,     (4.0+p)*L**2,  -6.0*L,     (2.0-p)*L**2],
            [-12.0,     -6.0*L,          12.0,     -6.0*L],
            [ 6.0*L,     (2.0-p)*L**2,  -6.0*L,     (4.0+p)*L**2]])

    def _mass_plane(self):
        """Consistent translational mass (includes shear effect through Phi)"""
        rho, A, L, p = self.material.rho, self.area, self.L, self.phi
        m11 = 13/35 + 7*p/10 + p**2/3
        m12 = (11/210 + 11*p/210 + p**2/24) * L
        m13 = 9/70 + 3*p/10 + p**2/6
        m14 = -(13/420 + 3*p/40 + p**2/24) * L
        m22 = (1/105 + p/60 + p**2/120) * L**2
        m24 = -(1/140 + p/60 + p**2/120) * L**2
        M = numpy.array([
            [m11,  m12,  m13,  m14],
            [m12,  m22, -m14,  m24],
            [m13, -m14,  m11, -m12],
            [m14,  m24, -m12,  m22]])
        return rho * A * L / (1.0 + p)**2 * M

    def _rotary_plane(self):
        """Consistent rotary-inertia matrix (includes shear effect through Phi)"""
        rho, I, L, p = self.material.rho, self.inertia, self.L, self.phi
        a = (3.0 - 15.0*p) * L
        d = (4.0 + 5.0*p + 10.0*p**2) * L**2
        e = (-1.0 - 5.0*p + 5.0*p**2) * L**2
        R = numpy.array([
            [ 36.0,  a,   -36.0,  a],
            [ a,     d,   -a,     e],
            [-36.0, -a,    36.0, -a],
            [ a,     e,   -a,     d]])
        return rho * I / (30.0 * L * (1.0 + p)**2) * R

    # Assembly of the two bending planes into the 8x8 element matrices
    def _assemble(self, plane_matrix):
        T = self._T
        out = numpy.zeros((8, 8))
        out[numpy.ix_(self._IDX_A, self._IDX_A)] = plane_matrix
        out[numpy.ix_(self._IDX_B, self._IDX_B)] = T @ plane_matrix @ T
        return out

    # Public element matrices
    def stiffness_matrix(self):
        """8x8 element stiffness matrix K"""
        return self._assemble(self._stiffness_plane())

    def mass_matrix(self):
        """8x8 element mass matrix M = M_translational (+ M_rotary if rotary_inertia)"""
        M = self._assemble(self._mass_plane())
        if self.rotary_inertia:
            M += self._assemble(self._rotary_plane())
        return M

    def gyroscopic_matrix(self):
        """
        8x8 gyroscopic matrix G (skew-symmetric), to be used as  C_gyro = Omega * G.
        Zero if gyroscopic=False.  Uses Ip = 2 Id, so G_AB = 2 R T, G_BA = -2 T R.
        """
        G = numpy.zeros((8, 8))
        if not self.gyroscopic:
            return G
        R, T = self._rotary_plane(), self._T
        G[numpy.ix_(self._IDX_A, self._IDX_B)] = 2.0 * R @ T
        G[numpy.ix_(self._IDX_B, self._IDX_A)] = -2.0 * T @ R
        return G

class Coupling:
    def __init__(self, shaft1:ShaftElement, shaft2:ShaftElement):
        self.shaft1 = shaft1
        self.shaft2 = shaft2

class Bearing:
    def __init__(self, damping, stiffness):
        self.C = damping
        self.k = stiffness

class Inertia:
    def __init__(self, inertia_matrix):
        self.I = inertia_matrix

class Rotor:
    def __init__(self, inertias, elements, bearings):
        self.inertias = inertias
        self.shaft_elements = elements
        self.bearings = bearings