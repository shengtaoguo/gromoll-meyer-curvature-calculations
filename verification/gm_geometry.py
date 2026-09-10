"""Reconstruct the exact A-family curvature model from its source geometry.

Public interface
----------------
    model = build_model(progress=lambda phase: None)

The returned dictionary has format ``GM-A-affine-mass-dual-v1``.  Its
20-dimensional matrix is A0+m*A1, its drive is s*(d0+m*d1), where
s^2=1-m, and its base is quadratic in m with rational and sqrt(2) parts.
The normalization is 15/2 times the sectional-curvature Taylor coefficient
before adding the feasible dual gain.  A rational polynomial pool represents
the angular coefficients.  No project files, saved formulas, or numerical
samples are read; this module does not write files.

The calculation retains metric-dependent horizontality and the entire
point-and-plane first-normal covector.  It verifies the explicit source
dual, its five homogeneous directions, the bracket-image projection, and
the complete symbolic mass identities.  It does not invert a normal
curvature Hessian.  Positivity of a chosen dual witness is a separate task
for the companion verifier.

Source provenance: the formulas were consolidated from the September 2026
research implementations listed in SOURCE_PROVENANCE below.  Old command
line wrappers, file bindings, numerical probes, B-family specializations,
and the eliminated normal-system solver are intentionally absent.  Each
remaining formula is ordinary Python source, rather than an encoded source
archive.  Dependencies: numpy, sympy, python-flint.
"""

import threading

import flint
import numpy as np
from sympy.polys.domains import QQ
from sympy.polys.fields import FracElement, field
from sympy.polys.rings import PolyElement


SOURCE_PROVENANCE = (
    "derive_a_mass_polynomial_model.py",
    "supplied_a_reverse_gradient.py",
    "supplied_a_explicit_dual_source.py",
    "probe_a_four_variable_reduction.py",
    "supplied_a_variational_geometry.py",
    "prove_supplied_a_first_null.py",
    "supplied_b_symbolic_geometry.py",
    "flint_exact_polynomial_backend.py",
)

_ARITHMETIC_LOCK = threading.RLock()


class _PolynomialBackend:
    """Temporarily accelerate exact field arithmetic, restoring all methods.

    SymPy owns the rational-function field.  Only polynomial multiplication
    and cancellation use FLINT, with exact cross checks against the original
    operations before the geometric calculation.  The NumPy dispatch guards
    prevent a rational-field scalar from trying to coerce an entire array.
    """

    def __init__(self, rational_field):
        self.field = rational_field
        self.contexts = {}
        self.counts = {"flint_products": 0, "flint_cancellations": 0}
        self.saved_fraction_methods = {}

    def context(self, ring):
        names = tuple(map(str, ring.symbols))
        if names not in self.contexts:
            self.contexts[names] = flint.fmpq_mpoly_ctx.get(names, "lex")
        return self.contexts[names]

    def native(self, poly):
        return self.context(poly.ring).from_dict(
            {mon: flint.fmpq(str(coefficient)) for mon, coefficient in poly.items()})

    @staticmethod
    def from_native(poly, ring):
        return ring.from_dict({mon: ring.domain.convert(c) for mon, c in poly.to_dict().items()})

    def multiply(self, first, second):
        if (isinstance(second, PolyElement) and first.ring == second.ring
                and first.ring.domain.is_QQ and first.ring.ngens >= 2
                and len(first)*len(second) > 256):
            self.counts["flint_products"] += 1
            return self.from_native(self.native(first)*self.native(second), first.ring)
        return self.original_multiply(first, second)

    def cancel(self, numerator, denominator):
        if (numerator.ring.domain.is_QQ and numerator.ring.ngens >= 2
                and len(numerator) > 1 and len(denominator) > 1):
            self.counts["flint_cancellations"] += 1
            nn, dd = self.native(numerator), self.native(denominator)
            common = nn.gcd(dd)
            nn, dd = nn/common, dd/common
            leading = dd.leading_coefficient()
            return (self.from_native(nn/leading, numerator.ring),
                    self.from_native(dd/leading, denominator.ring))
        return self.original_cancel(numerator, denominator)

    def restore(self):
        PolyElement.__mul__ = self.original_multiply
        PolyElement.cancel = self.original_cancel
        for name, method in self.saved_fraction_methods.items():
            setattr(FracElement, name, method)

    def __enter__(self):
        _ARITHMETIC_LOCK.acquire()
        self.original_multiply = PolyElement.__mul__
        self.original_cancel = PolyElement.cancel
        try:
            for name in ("__mul__", "__rmul__", "__add__", "__radd__",
                         "__sub__", "__rsub__", "__truediv__", "__rtruediv__"):
                method = getattr(FracElement, name)
                self.saved_fraction_methods[name] = method

                def guard(element, other, _original=method):
                    if isinstance(other, np.ndarray):
                        return NotImplemented
                    return _original(element, other)

                setattr(FracElement, name, guard)
            x, y, z = self.field.ring.gens
            cases = []
            for degree in range(1, 6):
                common = (x+y+z+1)**degree
                first = common*(x*x+2*y*z+3)
                second = common*(2*x*z+y*y+5)
                old_n, old_d = self.original_cancel(first, second)
                new_n, new_d = self.cancel(first, second)
                assert self.original_multiply(old_n, new_d) == self.original_multiply(new_n, old_d), "exact cancellation calibration"
                assert self.multiply(first, second) == self.original_multiply(first, second), "exact product calibration"
                cases.append((first, second, old_n, old_d))
            # Class attributes receive ordinary functions, not bound methods.
            backend = self
            PolyElement.__mul__ = lambda first, second: backend.multiply(first, second)
            PolyElement.cancel = lambda first, second: backend.cancel(first, second)
            for first, second, old_n, old_d in cases:
                actual = self.field.new(first, second)
                assert self.original_multiply(actual.numer, old_d) == self.original_multiply(old_n, actual.denom)
            a, b, c = self.field.gens
            assert (a+b)/(a-b)+(a-b)/(a+b)-2*(a*a+b*b)/(a*a-b*b) == 0
            assert (a+b+c)**3/(a+b+c)**2-(a+b+c) == 0
            self.calibration = {"exact_product_cases": 5, "exact_cancellation_cases": 5,
                                "rational_field_identities": 2, "flint_version": flint.__version__}
            return self
        except BaseException:
            self.restore()
            _ARITHMETIC_LOCK.release()
            raise

    def __exit__(self, exc_type, exc_value, traceback):
        self.restore()
        _ARITHMETIC_LOCK.release()
        return False


class _ExactEnvironment:
    """One field and its array constructors; no shared geometric globals."""

    def __init__(self):
        self.F, self.u, self.v, self.w = field("u,v,w", QQ)
        self.zero_scalar, self.one = self.F.zero, self.F.one
        self.algebra = None

    def Q(self, numerator, denominator=1):
        return self.F(QQ(numerator, denominator))

    def zeros(self, shape):
        return np.full(shape, self.zero_scalar, dtype=object)

    def eye(self, size):
        out = self.zeros((size, size))
        for i in range(size):
            out[i, i] = self.one
        return out

    @staticmethod
    def check_zero(value, label):
        assert not any(np.asarray(value, dtype=object).flat), label


def _left(q):
    """Real matrix for left multiplication by a quaternion."""
    a, b, c, d = q
    return np.array([[a, -b, -c, -d], [b, a, -d, c],
                     [c, d, a, -b], [d, -c, b, a]], dtype=object)


def _conjugate(q):
    return np.r_[q[0], -q[1:]]


def _quaternion_block(a, b, c, d):
    return np.block([[_left(a), _left(b)], [_left(c), _left(d)]])


def _cross(env, vector):
    a, b, c = vector
    zero = env.zero_scalar
    return np.array([[zero, -c, b], [c, zero, -a], [-b, a, zero]], dtype=object)


def _embed_k(env, vector):
    return np.r_[vector, env.zeros(4)]


def _embed_h(env, vector):
    return np.r_[vector, vector, env.zeros(4)]


class _Sp2Algebra:
    """The rational Lie frame: six diagonal and four off-diagonal coordinates."""

    def __init__(self, env):
        self.env = env
        e, zero = env.eye(4), env.zeros(4)
        self.basis = ([_quaternion_block(e[:, i], zero, zero, zero) for i in range(1, 4)]
                      + [_quaternion_block(zero, zero, zero, e[:, i]) for i in range(1, 4)]
                      + [_quaternion_block(zero, -_conjugate(e[:, i]), e[:, i], zero) for i in range(4)])
        self.q = np.diag([env.one]*6+[env.Q(2)]*4)
        adjoint = [np.column_stack([self.coordinates(A@B-B@A) for B in self.basis]) for A in self.basis]
        self.structure = [(k, i, j, adjoint[i][k, j])
                          for k in range(10) for i in range(10) for j in range(10)
                          if adjoint[i][k, j]]

    @staticmethod
    def coordinates(matrix):
        return np.array([matrix[1, 0], matrix[2, 0], matrix[3, 0],
                         matrix[5, 4], matrix[6, 4], matrix[7, 4],
                         matrix[4, 0], matrix[5, 0], matrix[6, 0], matrix[7, 0]], dtype=object)

    def matrix(self, vector):
        return sum((c*A for c, A in zip(vector, self.basis) if c), self.env.zeros((8, 8)))

    def bracket(self, first, second):
        out = self.env.zeros(10)
        for k, i, j, coefficient in self.structure:
            if first[i] and second[j]:
                out[k] += coefficient*first[i]*second[j]
        return out


def _product_bracket(env, first, second):
    """Bracket on Sp(2) x (Sp(1)xSp(1)) x Sp(1) x Sp(1)."""
    algebra = env.algebra
    out = env.zeros(22)
    out[:10] = algebra.bracket(first[:10], second[:10])
    out[10:16] = algebra.bracket(_embed_k(env, first[10:16]), _embed_k(env, second[10:16]))[:6]
    for offset in (16, 19):
        out[offset:offset+3] = env.Q(2)*_cross(env, first[offset:offset+3])@second[offset:offset+3]
    return out


def _phase(env, axis, parameter):
    denominator = env.one+env.Q(2)*parameter*parameter
    return np.r_[(env.one-env.Q(2)*parameter*parameter)/denominator,
                 env.Q(2)*parameter*axis/denominator]


def _rotation(env, quaternion):
    return np.column_stack([(_left(_left(quaternion)@np.r_[env.zero_scalar, env.eye(3)[:, i]])
                             @_conjugate(quaternion))[1:] for i in range(3)])


def _H_tilde(env, point):
    """The supplied H/sqrt(2), in the rational (not orthonormal) Lie frame."""
    Q, eye, zeros = env.Q, env.eye, env.zeros
    a, b = point[:4, 0], point[:4, 4]
    ra, rb, vb = a[0], b[0], b[1:]
    out = zeros((10, 10))
    out[3:6, 3:6] = Q(5)*ra*eye(3)
    out[:3, 7:] = -Q(4)*rb*eye(3)-Q(13, 15)*_cross(env, vb)
    out[3:6, 7:] = Q(107, 20)*rb*eye(3)-Q(4)*_cross(env, vb)
    out[3:6, 6] = Q(4)*vb
    out[7:, :3], out[7:, 3:6], out[6, 3:6] = out[:3, 7:].T, out[3:6, 7:].T, out[3:6, 6]
    return out


def _J_jets(env, points):
    """The original J=J_r+sqrt(2)J_i, with ordinary Taylor coefficients."""
    Q, zeros, eye, zero = env.Q, env.zeros, env.eye, env.zero_scalar
    a = np.array([point[:4, 0] for point in points], dtype=object)
    b = np.array([point[:4, 4] for point in points], dtype=object)
    mass = np.array([sum((a[i]@a[k-i] for i in range(k+1)), zero) for k in range(3)], dtype=object)
    real, imaginary = zeros((3, 10, 10)), zeros((3, 10, 10))
    for k in range(3):
        mass_squared = sum((mass[i]*mass[k-i] for i in range(k+1)), zero)
        real[k, :3, :3] = (Q(96)*mass[k]+Q(3)*mass_squared)*eye(3)
        for i in range(k+1):
            j = k-i
            real[k, 7:, 7:] += (-Q(56)*a[i, 0]*a[j, 0]*eye(3)
                                +Q(32)*np.outer(a[i, 1:], a[j, 1:])
                                -Q(98)*np.outer(b[i, 1:], b[j, 1:]))
            imaginary[k, 3:6, 7:] += Q(11)*b[i, 0]*_cross(env, a[j, 1:])
        imaginary[k, 7:, 3:6] = imaginary[k, 3:6, 7:].T
    return real, imaginary


def _flat_coefficients(env, point, X, Y, D, H, cu, cv):
    """Intrinsic J two-jets along the complete A flat generators."""
    algebra, Q, zero = env.algebra, env.Q, env.zero_scalar
    Xm, Ym, Dm, Hm = map(algebra.matrix, (X, Y, D, H))
    points = (point, point@(cu*Xm+cv*Ym),
              point@(cu*cu*(Xm@Xm)/Q(2)+cu*cv*(Xm@Ym)
                     +cv*cv*((Dm@Dm)/Q(2)+Dm@Hm+(Hm@Hm)/Q(2))))
    U = (X, -cv*algebra.bracket(H, X), cv*cv*algebra.bracket(H, algebra.bracket(H, X))/Q(2))
    V = (Y, -cv*algebra.bracket(H, D), cv*cv*algebra.bracket(H, algebra.bracket(H, D))/Q(2))

    def coefficient(first, jet, second):
        return sum((first[i]@jet[j]@second[2-i-j] for i in range(3) for j in range(3-i)), zero)

    return [tuple(coefficient(first, jet, second) for first, second in ((U, U), (U, V), (V, V)))
            for jet in _J_jets(env, points)]


class _ASource:
    """Actual product-submersion geometry at the complete rational A flag.

    The source metric has factors Q_G,Q_K,Q_H,(1/10)Q_L.  Its fixed
    horizontal commuting pair has Gram determinant 15.  The physical
    projection is (g,k,h) -> g k^-1 K1 h^-1, so its differential and all
    differentiated tensor terms retain the K and H contributions.
    """

    def __init__(self, env, pu, pv, pw, progress):
        self.env, self.a, self.progress = env, env.algebra, progress
        Q, zeros, eye, one, zero = env.Q, env.zeros, env.eye, env.one, env.zero_scalar
        y = np.array([one, one, zero], dtype=object)
        z = np.array([-one, zero, one], dtype=object)
        axis = y+z
        A0 = np.array([one, one, one, one], dtype=object)/Q(2)
        B0 = np.array([one, -one, one, one], dtype=object)/Q(2)
        qa, qb = _phase(env, y, pu), _phase(env, z, pv)
        r, s = (one-pw*pw)/(one+pw*pw), Q(2)*pw/(one+pw*pw)
        self.point = _quaternion_block(r*A0, s*B0, -s*A0, r*B0)
        self.K1 = _quaternion_block(qa, zeros(4), zeros(4), qb)
        self.physical = self.point@self.K1
        self.y, self.z, self.axis = y, z, axis
        env.check_zero(self.physical.T@self.physical-eye(8), "complete rational A family lies in Sp(2)")
        self.q = np.diag([one]*6+[Q(2)]*4+[one]*9+[Q(1, 10)]*3)
        self.qi = np.diag([one]*6+[Q(1, 2)]*4+[one]*9+[Q(10)]*3)
        U, V = zeros(22), zeros(22)
        U[6:10] = _left(_conjugate(B0))@A0
        V[:3], V[3:6] = y/Q(2), z/Q(2)
        V[10:13], V[13:16], V[16:19], V[19:22] = -y/Q(2), -z/Q(2), -axis/Q(2), Q(5)*axis
        self.U, self.V = U, V
        self.gram = (U@self.q@U)*(V@self.q@V)-(U@self.q@V)**2
        assert self.gram == Q(15), "the actual Cheeger-10 source Gram determinant"
        self.SA, self.SB = _rotation(env, qa), _rotation(env, qb)
        L, R = zeros((22, 15)), zeros((22, 15))
        for i in range(6):
            R[i, i] = R[10+i, i] = one
        L[10:13, 6:9], L[13:16, 6:9] = self.SA, self.SB
        for i in range(3):
            R[16+i, 6+i] = one
            L[3+i, 9+i] = R[19+i, 9+i] = one
            L[i, 12+i] = L[16+i, 12+i] = one
        for i in range(15):
            L[:10, i] = self.a.coordinates(self.point.T@self.a.matrix(L[:10, i])@self.point)
        self.L, self.R, self.vertical = L, R, L-R
        env.check_zero(self.vertical.T@self.q@U, "A source U is horizontal")
        env.check_zero(self.vertical.T@self.q@V, "A source V is horizontal")
        env.check_zero(_product_bracket(env, U, V), "the A source flat is commuting")
        self.AU = np.column_stack([_product_bracket(env, U, eye(22)[:, i]) for i in range(22)])
        self.AV = np.column_stack([_product_bracket(env, V, eye(22)[:, i]) for i in range(22)])
        self.T = np.column_stack([self.push(eye(22)[:, i]) for i in range(22)])
        env.check_zero(self.T@self.vertical[:, :9], "the physical differential kills the first two quotient fibres")
        self.hp = _H_tilde(env, self.physical)
        self.h = self.T.T@self.hp@self.T
        self.hU, self.hV = self.h@U, self.h@V
        self.dhU, self.dhV = self.first_metric(U), self.first_metric(V)
        self.Cuu = self.qi@(self.dhU@U-self.gradient_metric(U, U)/Q(2)-self.AU.T@self.hU)
        self.Cvv = self.qi@(self.dhV@V-self.gradient_metric(V, V)/Q(2)-self.AV.T@self.hV)
        self.Cuv = self.qi@(self.dhU@V+self.dhV@U-self.gradient_metric(U, V)
                           -self.AV.T@self.hU-self.AU.T@self.hV)/Q(2)
        self.fixed_group_second_coefficient = self.Cuv@self.q@self.Cuv-self.Cuu@self.q@self.Cvv
        self.lower_covector = Q(3, 2)*(self.dhU@V-self.dhV@U)+(self.AV.T@self.hU-self.AU.T@self.hV)/Q(2)
        self.omega_h = (-(self.AU@L).T@self.hV+(self.AV@L).T@self.hU
                        +self.vertical.T@(self.dhU@V-self.dhV@U))
        env.check_zero(self.vertical[:, :9].T@self.h, "H is basic for the first two quotient fibres")
        progress("A-connection-and-metric-ONeill-variation-ready")

    def adjK(self, vector):
        return self.a.coordinates(self.K1.T@self.a.matrix(vector)@self.K1)

    def push(self, vector):
        return self.adjK(vector[:10]-_embed_k(self.env, vector[10:16]))-_embed_h(self.env, vector[16:19])

    def dpush(self, direction, vector):
        env = self.env
        return (self.a.bracket(_embed_h(env, direction[16:19]), self.push(vector))
                +self.adjK(self.a.bracket(_embed_k(env, direction[10:16]),
                                          vector[:10]-_embed_k(env, vector[10:16]))))

    def first_metric(self, direction):
        env = self.env
        dT = np.column_stack([self.dpush(direction, env.eye(22)[:, i]) for i in range(22)])
        hp1 = _H_tilde(env, self.physical@self.a.matrix(self.push(direction)))
        return self.T.T@hp1@self.T+dT.T@self.hp@self.T+self.T.T@self.hp@dT

    def gradient_metric(self, first, second):
        env = self.env
        pushed_first, pushed_second = self.push(first), self.push(second)
        derivative = np.array([pushed_first@_H_tilde(env, self.physical@E)@pushed_second
                               for E in self.a.basis], dtype=object)
        DA = np.column_stack([self.dpush(env.eye(22)[:, i], first) for i in range(22)])
        DB = np.column_stack([self.dpush(env.eye(22)[:, i], second) for i in range(22)])
        return self.T.T@derivative+DA.T@self.hp@pushed_second+DB.T@self.hp@pushed_first

    def J_variation(self):
        env, Q = self.env, self.env.Q
        X, Y = self.push(self.U), self.push(self.V)
        D, H = env.zeros(10), env.zeros(10)
        D[:3], D[3:6] = self.y, self.z
        H[:3] = H[3:6] = self.axis/Q(2)
        values = [_flat_coefficients(env, self.physical, X, Y, D, H, cu, cv)
                  for cu, cv in ((env.one, env.one), (env.one, env.zero_scalar), (env.zero_scalar, env.one))]
        return [(values[0][i][1]-values[1][i][1]-values[2][i][1]-values[1][i][2]-values[2][i][0])/Q(15)
                for i in range(2)]


class _Node:
    """Matrix reverse differentiation with exact field-valued entries."""

    def __init__(self, graph, value, parents=(), variable=False):
        self.graph, self.value = graph, np.asarray(value, dtype=object)
        assert self.value.ndim == 2
        self.parents, self.active, self.gradient = parents, bool(variable or parents), None
        if self.active:
            graph.nodes.append(self)

    def __add__(self, other):
        assert self.value.shape == other.value.shape
        parents = []
        if self.active:
            parents.append((self, lambda derivative: derivative))
        if other.active:
            parents.append((other, lambda derivative: derivative))
        return _Node(self.graph, self.value+other.value, parents)

    def __sub__(self, other):
        return self+other.scaled(-self.graph.env.one)

    def scaled(self, coefficient):
        parents = [(self, lambda derivative: coefficient*derivative)] if self.active else []
        return _Node(self.graph, coefficient*self.value, parents)

    def __matmul__(self, other):
        parents = []
        if self.active:
            parents.append((self, lambda derivative: derivative@other.value.T))
        if other.active:
            parents.append((other, lambda derivative: self.value.T@derivative))
        return _Node(self.graph, self.value@other.value, parents)

    def transpose(self):
        parents = [(self, lambda derivative: derivative.T)] if self.active else []
        return _Node(self.graph, self.value.T, parents)


class _Graph:
    def __init__(self, env, progress):
        self.env, self.progress, self.nodes = env, progress, []

    def constant(self, value):
        return _Node(self, value)

    def variable(self, value):
        return _Node(self, value, variable=True)

    def linear(self, source, forward, adjoint):
        return _Node(self, forward(source.value), [(source, adjoint)] if source.active else [])

    def backward(self, output):
        assert output.value.shape == (1, 1)
        output.gradient = np.array([[self.env.one]], dtype=object)
        for index, node in enumerate(reversed(self.nodes)):
            if node.gradient is not None and any(node.gradient.flat):
                for parent, pullback in node.parents:
                    contribution = pullback(node.gradient)
                    parent.gradient = contribution.copy() if parent.gradient is None else parent.gradient+contribution
            if index % 50 == 0:
                self.progress(f"reverse-principal-node-{index+1}-of-{len(self.nodes)}")


def _principal_gradient(calc):
    """Differentiate the complete first-curvature symbol in all 54 variables.

    The three phase parameters are constants for this differentiation.
    The independent variables are the source point and both source vectors;
    the point gradient is subsequently contracted with all ten Sp(2) fields.
    """
    env, algebra = calc.env, calc.a
    Q, zeros, eye, zero, one = env.Q, env.zeros, env.eye, env.zero_scalar, env.one
    graph = _Graph(env, calc.progress)
    point = graph.variable(calc.point)
    U, V = graph.variable(calc.U.reshape((22, 1))), graph.variable(calc.V.reshape((22, 1)))
    K1, T = graph.constant(calc.K1), graph.constant(calc.T)
    difference, PK, PH, PG = (zeros((10, 22)) for _ in range(4))
    difference[:, :10], difference[:6, 10:16] = eye(10), -eye(6)
    PK[:6, 10:16] = eye(6)
    PH[:3, 16:19] = PH[3:6, 16:19] = eye(3)
    PG[:, :10] = eye(10)
    PG, PK, PH, difference = map(graph.constant, (PG, PK, PH, difference))
    adjK = graph.constant(np.column_stack([calc.adjK(eye(10)[:, i]) for i in range(10)]))
    positions = ((1, 0), (2, 0), (3, 0), (5, 4), (6, 4), (7, 4), (4, 0), (5, 0), (6, 0), (7, 0))

    def matrix(source):
        return graph.linear(source, lambda value: algebra.matrix(value[:, 0]),
                            lambda derivative: np.array([[sum((derivative*B).flat, zero)] for B in algebra.basis], dtype=object))

    def coordinates(source):
        def adjoint(derivative):
            out = zeros((8, 8))
            for i, (j, k) in enumerate(positions):
                out[j, k] += derivative[i, 0]
            return out
        return graph.linear(source, lambda value: algebra.coordinates(value).reshape((10, 1)), adjoint)

    def bracket(first, second):
        first_matrix, second_matrix = matrix(first), matrix(second)
        return coordinates(first_matrix@second_matrix-second_matrix@first_matrix)

    slots = ((0, 0), (0, 4), (1, 4), (2, 4), (3, 4))
    coefficients = []
    for i, j in slots:
        basis = zeros((8, 8)); basis[i, j] = one
        coefficients.append(_H_tilde(env, basis))

    def hmap(source):
        def adjoint(derivative):
            out = zeros((8, 8))
            for (i, j), coefficient in zip(slots, coefficients):
                out[i, j] = sum((derivative*coefficient).flat, zero)
            return out
        return graph.linear(source, lambda value: _H_tilde(env, value), adjoint)

    def inner(first, metric, second):
        return first.transpose()@metric@second

    physical = point@K1
    h0 = hmap(physical)
    pushed_cache, flow_cache = {}, {}

    def pushed(direction, vector):
        key = id(direction), id(vector)
        if key in pushed_cache:
            return pushed_cache[key]
        value, HW, KW, DD = T@vector, PH@direction, PK@direction, difference@vector
        kfirst = adjK@bracket(KW, DD)
        first = bracket(HW, value)+kfirst
        second = (bracket(HW, bracket(HW, value))+bracket(HW, kfirst).scaled(Q(2))
                  +adjK@bracket(KW, bracket(KW, DD)))
        pushed_cache[key] = value, first, second
        return value, first, second

    def flow(direction):
        if id(direction) in flow_cache:
            return flow_cache[id(direction)]
        WG, WK, WH = matrix(PG@direction), matrix(PK@direction), matrix(PH@direction)
        first = point@(WG-WK)@K1-physical@WH
        second = (point@(WG@WG-(WG@WK).scaled(Q(2))+WK@WK)@K1
                  -(point@(WG-WK)@K1@WH).scaled(Q(2))+physical@WH@WH)
        flow_cache[id(direction)] = hmap(first), hmap(second)
        return flow_cache[id(direction)]

    def second(direction, first_vector, second_vector):
        A0, A1, A2 = pushed(direction, first_vector)
        B0, B1, B2 = pushed(direction, second_vector)
        h1, h2 = flow(direction)
        return (inner(A0, h2, B0)+inner(A1, h1, B0).scaled(Q(2))+inner(A0, h1, B1).scaled(Q(2))
                +inner(A2, h0, B0)+inner(A1, h0, B1).scaled(Q(2))+inner(A0, h0, B2))

    output = (second(U+V, U, V)-second(U, U, V)-second(V, U, V)
              -second(U, V, V)-second(V, U, U)).scaled(Q(1, 2))
    assert output.value[0, 0] == 0, "global A first-nullity in the full source expression"
    calc.progress(f"principal-graph-built-{len(graph.nodes)}-active-nodes")
    graph.backward(output)
    calc.bg = np.array([sum((point.gradient*(calc.point@E)).flat, zero) for E in algebra.basis], dtype=object)
    calc.bu, calc.bv = U.gradient[:, 0], V.gradient[:, 0]
    calc.progress("complete-point-and-plane-principal-covector-ready")


def _raw_constraints(calc):
    """E delta+f=0 and beta=B delta, with all 30 horizontal equations."""
    env = calc.env
    E = env.zeros((30, 54))
    for i in range(10):
        dvertical = env.zeros((22, 15))
        dvertical[:10] = np.column_stack([-calc.a.bracket(env.eye(10)[:, i], calc.L[:10, j]) for j in range(15)])
        E[:15, i] = dvertical.T@calc.q@calc.U
        E[15:, i] = dvertical.T@calc.q@calc.V
    E[:15, 10:32] = calc.vertical.T@calc.q
    E[15:, 32:] = calc.vertical.T@calc.q
    B = env.zeros((22, 54))
    B[:, 10:32], B[:, 32:] = -calc.AV, calc.AU
    kappa = np.r_[calc.bg, calc.bu, calc.bv]+B.T@calc.lower_covector
    forcing = np.r_[calc.vertical.T@calc.hU, calc.vertical.T@calc.hV]
    return E, B, kappa, forcing


def _explicit_source_dual(calc, kappa):
    """A particular dual and five homogeneous directions with constant divisors."""
    env, algebra = calc.env, calc.a
    Q, zeros, eye, one, zero = env.Q, env.zeros, env.eye, env.one, env.zero_scalar
    y, z, axis = calc.y, calc.z, calc.axis
    normal = _cross(env, y)@z
    A0 = np.array([one]*4, dtype=object)/Q(2)
    B0 = np.array([one, -one, one, one], dtype=object)/Q(2)
    T0 = _quaternion_block(A0, zeros(4), zeros(4), B0)

    def Ad(matrix, vector):
        return algebra.coordinates(matrix@algebra.matrix(vector)@matrix.T)

    RA0, RB0 = _rotation(env, A0), _rotation(env, B0)
    SA, SB = calc.SA, calc.SB
    kU, kV = calc.qi@kappa[10:32], calc.qi@kappa[32:]
    leftU, leftV = Ad(T0, kU[:10]), Ad(T0, kV[:10])
    gamma = Ad(calc.point, calc.qi[:10, :10]@kappa[:10])
    r, s = Q(2)*calc.point[0, 0], Q(2)*calc.point[0, 4]
    RA = RA0@SA
    qperp = -(_cross(env, axis)@gamma[:3])/Q(2)
    first = leftV[:3]-RA@kV[16:19]-RA0@kV[10:13]
    w0 = (first+(RA-r*r*eye(3))@qperp)/Q(2)
    zA = _cross(env, axis)@RA@axis
    env.check_zero(zA@zA-Q(3), "the explicit source divisor is the constant three")
    tau = Q(2, 3)*(zA@(leftU[7:]-r*s*gamma[7:]-_cross(env, axis)@w0))
    qV = qperp+tau*axis
    lamV = zeros(15)
    lamV[12:], lamV[6:9] = qV, qV-kV[16:19]
    lamV[:3], lamV[3:6] = SA@lamV[6:9]-kV[10:13], SB@lamV[6:9]-kV[13:16]
    difference = leftV[7:]-r*s*qV
    zg = zeros(10)
    zg[:3], zg[3:6] = -difference/Q(2), difference/Q(2)
    zg[7:] = (first+(RA-r*r*eye(3))@qV)/Q(2)
    Z = zeros(22)
    Z[:10] = Ad(T0.T, zg)
    lamU = zeros(15)
    lamU[12:] = gamma[7:]
    lamU[:3] = RA0.T@(_cross(env, axis)@zg[:3]+r*r*gamma[7:]-leftU[:3])
    lamU[3:6] = RB0.T@(_cross(env, axis)@zg[3:6]+s*s*gamma[7:]-leftU[3:6])
    b1, b2 = y@(kU[10:13]+lamU[:3]), z@(kU[13:16]+lamU[3:6])
    lamU[6:9] = ((Q(2)*b1+b2)*y+(b1+Q(2)*b2)*z)/Q(3)
    Z[10:13] = _cross(env, y)@(kU[10:13]+lamU[:3]-SA@lamU[6:9])/Q(2)
    Z[13:16] = _cross(env, z)@(kU[13:16]+lamU[3:6]-SB@lamU[6:9])/Q(2)
    Z[16:19] = _cross(env, axis)@(kU[16:19]+lamU[6:9]-gamma[7:])/Q(2)
    lam = np.r_[lamU, lamV]
    kernel_Z, kernel_lambda = [], []
    # Coordinates: e_n,e_(y-z),t_n,t_(y-z),nu, where e=axis x zfree+t.
    for e, t, nu in ((normal, zeros(3), zero), (y-z, zeros(3), zero),
                     (zeros(3), normal, zero), (zeros(3), y-z, zero),
                     (zeros(3), zeros(3), one)):
        zfree = -(_cross(env, axis)@(e-t))/Q(2)
        q1, q2, h0 = RA0.T@e, RB0.T@e, nu*normal
        kk, gg = zeros(22), zeros(10)
        gg[:3] = gg[3:6] = zfree
        kk[:10] = Ad(T0.T, gg)
        kk[10:13] = _cross(env, y)@(q1-SA@h0)/Q(2)
        kk[13:16] = _cross(env, z)@(q2-SB@h0)/Q(2)
        kk[16:19] = _cross(env, axis)@(h0-t)/Q(2)
        kk[19:] = -(_cross(env, axis)@t)/Q(20)
        kernel_Z.append(kk)
        kernel_lambda.append(np.r_[q1, q2, h0, t, t, zeros(15)])
    return Z, lam, np.column_stack(kernel_Z), np.column_stack(kernel_lambda)


def _twenty_variable_model(calc, E, B, kappa):
    """Expand a feasible source dual, then use a constant coordinate change.

    P is the Q-orthogonal projection onto the bracket image.  For variables
    (z,eta), use Z+KZ*z-P*Q^-1*C^T*eta and lambda+KL*z in the completed
    square formula.  The returned matrix and drive are its exact quadratic
    penalty and linear gain, after the displayed constant 20-by-20 change.
    """
    env = calc.env
    Q, zeros, eye, one, zero = env.Q, env.zeros, env.eye, env.one, env.zero_scalar
    Z, lam, KZ, KL = _explicit_source_dual(calc, kappa)
    env.check_zero(B.T@calc.q@Z+E.T@lam-kappa, "all 54 inhomogeneous source-dual equations")
    calc.progress("all-54-source-dual-equations-verified")
    env.check_zero(B.T@calc.q@KZ+E.T@KL, "all five homogeneous source-dual vectors")
    calc.progress("all-five-homogeneous-source-duals-verified")
    AU, AV = calc.AU, calc.AV
    PU = -AU[:10, :10]@AU[:10, :10]/Q(4)
    PV = -AV[:10, :10]@AV[:10, :10]/Q(2)
    P = zeros((22, 22))
    P[:10, :10] = PU+PV-PU@PV
    for offset in (10, 13, 16, 19):
        axis = calc.V[offset:offset+3]
        P[offset:offset+3, offset:offset+3] = eye(3)-np.outer(axis, axis)/(axis@axis)
    for value, label in ((P@P-P, "bracket-image projection is idempotent"),
                         (P.T@calc.q-calc.q@P, "bracket-image projection is Q self-adjoint"),
                         (P@B-B, "the projection fixes the full bracket map"),
                         (P@Z-Z, "the particular dual lies in the bracket image"),
                         (P@KZ-KZ, "the homogeneous duals lie in the bracket image")):
        env.check_zero(value, label)
    C = (calc.L+calc.R).T@calc.q
    G = calc.vertical.T@calc.q@calc.vertical
    Kg, coupling = KZ.T@calc.q@KZ, KZ.T@C.T
    Aeta = C@P@calc.qi@C.T+G/Q(3)
    joint = np.block([[Kg, -coupling], [-coupling.T, Aeta]])
    env.check_zero(joint-joint.T, "complete dual quadratic matrix is symmetric")
    A0 = np.array([one]*4, dtype=object)/Q(2)
    B0 = np.array([one, -one, one, one], dtype=object)/Q(2)
    RA0, RB0 = _rotation(env, A0), _rotation(env, B0)
    axis, normal, perpendicular = calc.axis, _cross(env, calc.y)@calc.z, calc.y-calc.z
    transform = zeros((20, 20))

    def group_column(column, direction, part):
        offset = (0, 3, 9, 12)[part]
        value = RA0.T@direction if part == 0 else RB0.T@direction if part == 1 else direction
        transform[5+offset:5+offset+3, column] = value

    for j in range(4):
        group_column(j, axis, j)
    for direction, begin, ep, tp in ((normal, 4, 1, 3), (perpendicular, 10, 0, 2)):
        for j in range(4):
            group_column(begin+j, direction, j)
        transform[ep, begin+4] = -Q(1, 2) if begin == 4 else one
        transform[tp, begin+5] = -Q(1, 2) if begin == 4 else one
    transform[4, 16] = one
    transform[11:14, 17:20] = eye(3)
    rotated = transform.T@joint@transform
    env.check_zero(rotated[:4, 4:16], "parallel bulk block decouples")
    env.check_zero(rotated[4:10, 10:16], "the two transverse bulk blocks decouple")
    env.check_zero(rotated[4:10, 4:10]/Q(3)-rotated[10:16, 10:16]/Q(6), "identical transverse scalar blocks")
    calc.progress("exact-20-dimensional-dual-model-ready")
    return Z, lam, KZ, KL, C, G, transform, rotated


def _assemble(env, parameter, label, progress):
    def report(phase):
        progress(label+"-"+phase)

    calc = _ASource(env, env.u, env.v, parameter, report)
    _principal_gradient(calc)
    E, B, kappa, forcing = _raw_constraints(calc)
    Z, lam, KZ, KL, C, G, transform, matrix = _twenty_variable_model(calc, E, B, kappa)
    drive = transform.T@np.r_[-KZ.T@calc.q@Z-KL.T@forcing/env.Q(2),
                              C@Z-calc.omega_h/env.Q(2)]
    base = calc.fixed_group_second_coefficient-lam@forcing-Z@calc.q@Z
    Jr, Ji = calc.J_variation()
    report("source-and-intrinsic-J-assembled")
    return {"matrix": matrix, "drive": drive,
            "baseR": base+env.Q(15, 2)*Jr, "baseI": env.Q(15, 2)*Ji}


def _encode_model(env, A0, A1, d0, d1, base, checks):
    """Return rational coefficient data in memory; no file operations."""
    pool, indices = [], {}

    def poly(value):
        assert all(mon[2] == 0 for mon in value), "every model coefficient depends only on the two phases"
        key = frozenset(value.items())
        if key not in indices:
            indices[key] = len(pool)
            pool.append([[int(mon[0]), int(mon[1]), str(c)] for mon, c in sorted(value.items())])
        return indices[key]

    def scalar(value):
        return [poly(value.numer), poly(value.denom)]

    def array(values):
        data = np.asarray(values, dtype=object)
        return scalar(data.item()) if data.ndim == 0 else [array(row) for row in data]

    return {"format": "GM-A-affine-mass-dual-v1", "angle_symbols": ["u", "v"],
            "mass_interval": [0, 1], "matrix": "A0+mass*A1",
            "drive": "sqrt(1-mass)*(d0+mass*d1)",
            "drive_representative": "s>=0",
            "signed_drive": "s*(d0+mass*d1)",
            "row_factor_relation": "s^2=1-mass",
            "base_normalization": "15/2 times the sectional Taylor coefficient before adding dual gain",
            "retained_order": ["nu", "eta_h_i", "eta_h_j", "eta_h_k"],
            "source_provenance": list(SOURCE_PROVENANCE), "exact_checks": checks,
            "polynomials": pool, "A0": array(A0), "A1": array(A1),
            "d0": array(d0), "d1": array(d1),
            "baseR": array(base["baseR"]), "baseI": array(base["baseI"])}


def build_model(progress=lambda phase: None):
    """Rebuild and verify the whole A mass model over QQ(u,v,w).

    The three row values construct candidate coefficients only.  The last
    assembly uses an independent symbolic row variable and verifies all
    entries on the entire family, so interpolation samples are not evidence
    for the mass identities.  Both row poles remain explicit inputs.
    """
    if not __debug__:
        raise RuntimeError("exact verification requires assertions; do not run with -O")
    assert callable(progress)
    env = _ExactEnvironment()
    with _PolynomialBackend(env.F) as backend:
        env.algebra = _Sp2Algebra(env)
        Q, one, zero = env.Q, env.one, env.zero_scalar
        progress("exact-polynomial-backend-calibrated")
        mass_zero = _assemble(env, one, "mass-zero", progress)
        progress("angle-only-mass-zero-source-complete")
        mass_one = _assemble(env, zero, "mass-one", progress)
        env.check_zero(mass_one["drive"], "the full dual drive vanishes at the b-zero row pole")
        progress("angle-only-mass-one-source-and-drive-vanishing-complete")
        middle = _assemble(env, Q(1, 2), "mass-nine-over-twenty-five", progress)
        middle_mass = Q(9, 25)
        d0 = mass_zero["drive"]
        d1 = (middle["drive"]/Q(4, 5)-d0)/middle_mass
        A0, A1 = mass_zero["matrix"], mass_one["matrix"]-mass_zero["matrix"]
        base = {}
        for key in ("baseR", "baseI"):
            c0 = mass_zero[key]
            c2 = (middle[key]-(one-middle_mass)*c0-middle_mass*mass_one[key])/(middle_mass*(middle_mass-one))
            c1 = mass_one[key]-c0-c2
            base[key] = [c0, c1, c2]
        progress("candidate-affine-and-quadratic-mass-coefficients-constructed")
        full = _assemble(env, env.w, "complete-symbolic-row-family", progress)
        row = (one-env.w*env.w)/(one+env.w*env.w)
        s = Q(2)*env.w/(one+env.w*env.w)
        mass = row*row
        assert mass+s*s == one
        for i in range(20):
            env.check_zero(full["matrix"][i]-A0[i]-mass*A1[i], "full global affine mass matrix row "+str(i))
            env.check_zero(full["drive"][i]-s*(d0[i]+mass*d1[i]), "full global row-factor drive entry "+str(i))
            progress("complete-mass-matrix-and-drive-row-"+str(i+1)+"-verified")
        for key in ("baseR", "baseI"):
            env.check_zero(full[key]-base[key][0]-mass*base[key][1]-mass*mass*base[key][2],
                           "full global quadratic mass base "+key)
            progress("complete-quadratic-mass-"+key+"-verified")
        for entry in A1.flat:
            assert all(sum(mon) == 0 for mon in entry.numer) and all(sum(mon) == 0 for mon in entry.denom), "the mass matrix slope is constant"
        checks = {"global_affine_matrix_exact": True, "global_row_factor_drive_exact": True,
                  "global_quadratic_base_exact": True, "A_b_zero_row_pole_drive_exact": True,
                  "mass_matrix_slope_constant": True, "global_first_null_exact": True,
                  "full_source_dual_and_five_kernel_identities_exact": True,
                  "bracket_image_projection_identities_exact": True,
                  "metric_dependent_horizontality_retained": True,
                  "gram_determinant": "15", "backend_calibration": backend.calibration}
        model = _encode_model(env, A0, A1, d0, d1, base, checks)
        progress("geometry-to-20-dimensional-mass-model-complete")
        return model
