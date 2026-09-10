"""Run the GM construction's final exact algebraic checks.

Usage: python verify_gm.py
Keep gm_geometry.py beside this file. No other project files or saved
polynomial/model/certificate data are read. Dependencies: numpy, sympy,
python-flint. Optional --output and --status write machine-readable results;
the default run writes no files and prints a short progress report.

gm_geometry.py rebuilds the source model from the metric formulas. This
independent checker uses native FLINT polynomials and Fraction arithmetic
to reconstruct the witness and prove its five angular signs. The finite
partition below is a proposed cover: every bound and both children of
every split are checked afresh. No saved numerator or sign is trusted.

The accompanying mathematical argument supplies the quotient interpretation,
complete zero-stratum coverage, smoothness, and compactness completion.
This script checks the algebraic part; it is not a formalization of those
geometric lemmas.
"""

import argparse
import hashlib
import json
import math
import platform
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path

sys.dont_write_bytecode = True

try:
    import flint
except ImportError as error:
    raise SystemExit("Missing dependency. Run: python -m pip install numpy sympy python-flint") from error


CTX = flint.fmpq_mpoly_ctx.get(("u", "v"), "lex")
ZERO, ONE = Q(0), Q(1)

# (split axis, coefficient indices proved positive at this node).
# -1 denotes a leaf. Nodes are in preorder; the verifier consumes both
# children and rejects uncovered functions or unused nodes.
PARTITION = [
    (1, ()), (0, ()), (0, (0, 1, 2, 3)),
    (1, ()), (-1, (4,)), (-1, (4,)),
    (1, ()), (-1, (4,)), (-1, (4,)),
    (0, (0, 1, 2, 3)),
    (1, ()), (-1, (4,)), (-1, (4,)),
    (1, ()), (-1, (4,)), (-1, (4,)),
    (0, ()), (1, (2, 3)), (-1, (0, 1, 4)), (-1, (0, 1, 4)),
    (1, (3,)), (0, (0, 4)), (1, ()),
    (-1, (1, 2)), (-1, (1, 2)), (-1, (1, 2)),
    (0, (0, 1, 2)), (-1, (4,)), (-1, (4,)),
]

WITNESS = {
    "reference_mass": "1/3", "radial_bernstein_degree": 4,
    "filter_coefficients": ["17/8", "-21/16", "1/4"],
    "sqrt2_rational_comparison": "7/5", "rational_margin": "1/8",
    "angular_box": [["-3/4", "3/4"], ["-3/4", "3/4"]],
}


def fp(value):
    value = Q(value)
    return flint.fmpq(value.numerator, value.denominator)


def polynomial(terms):
    return CTX.from_dict({tuple(mon): fp(c) for mon, c in terms.items() if c})


def constant(value):
    return polynomial({(0, 0): Q(value)})


DU = polynomial({(2, 0): ONE, (0, 0): Q(1, 2)})
DV = polynomial({(0, 2): ONE, (0, 0): Q(1, 2)})


@lru_cache(maxsize=None)
def circle_power(a, b):
    assert a >= 0 and b >= 0
    return DU**a * DV**b


def coefficient_dict(poly):
    return {tuple(map(int, mon)): Q(str(c)) for mon, c in poly.to_dict().items()}


def polynomial_degrees(poly):
    monomials = poly.to_dict()
    return tuple(max((int(mon[j]) for mon in monomials), default=0) for j in (0, 1))


class RatCircle:
    """Immutable numerator / DU^a DV^b; arithmetic never performs a gcd."""

    __slots__ = ("n", "a", "b")

    def __init__(self, numerator, a=0, b=0):
        assert isinstance(a, int) and isinstance(b, int) and min(a, b) >= 0
        self.n = numerator
        self.a, self.b = (a, b) if numerator else (0, 0)

    @staticmethod
    def scalar(value):
        return RatCircle(constant(Q(value)))

    @staticmethod
    def coerce(value):
        return value if isinstance(value, RatCircle) else RatCircle.scalar(value)

    def __bool__(self):
        return bool(self.n)

    def lifted(self, a, b):
        assert a >= self.a and b >= self.b
        return self.n * circle_power(a-self.a, b-self.b)

    def __add__(self, other):
        other = self.coerce(other)
        if not self:
            return other
        if not other:
            return self
        a, b = max(self.a, other.a), max(self.b, other.b)
        return RatCircle(self.lifted(a, b)+other.lifted(a, b), a, b)

    __radd__ = __add__

    def __neg__(self):
        return RatCircle(-self.n, self.a, self.b)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        return RatCircle(self.n*other.n, self.a+other.a, self.b+other.b)

    __rmul__ = __mul__

    def __truediv__(self, other):
        # General rational-function division is deliberately unavailable.
        assert not isinstance(other, RatCircle)
        other = Q(other)
        assert other
        return RatCircle(self.n*fp(1/other), self.a, self.b)

    def __eq__(self, other):
        other = self.coerce(other)
        a, b = max(self.a, other.a), max(self.b, other.b)
        return self.lifted(a, b) == other.lifted(a, b)

    def cross_equal(self, other):
        other = self.coerce(other)
        return (self.n*circle_power(other.a, other.b)
                == other.n*circle_power(self.a, self.b))

    def as_fraction(self):
        value = Q(str(self.n.to_dict().get((0, 0), 0))) * Q(2)**(self.a+self.b)
        assert self.n == constant(value)*circle_power(self.a, self.b), "the claimed constant must be constant at every angle"
        return value

    def at(self, u, v):
        numerator = sum((c*u**i*v**j for (i, j), c in coefficient_dict(self.n).items()), ZERO)
        return numerator/(u*u+Q(1, 2))**self.a/(v*v+Q(1, 2))**self.b


def circle_denominator(denominator):
    assert denominator, "zero model denominator"
    du, dv = polynomial_degrees(denominator)
    assert du % 2 == 0 and dv % 2 == 0, "circle denominator has even bidegree"
    lead = Q(str(denominator.to_dict().get((du, dv), 0)))
    assert lead and denominator == constant(lead)*circle_power(du//2, dv//2), "complete exact recognition of the positive circle factors"
    return lead, du//2, dv//2


def decode_model(source):
    assert source["format"] == "GM-A-affine-mass-dual-v1"
    assert source["angle_symbols"] == ["u", "v"]
    assert source["mass_interval"] == [0, 1]
    assert source["matrix"] == "A0+mass*A1"
    assert source["drive"] == "sqrt(1-mass)*(d0+mass*d1)"
    assert source["drive_representative"] == "s>=0"
    assert source["signed_drive"] == "s*(d0+mass*d1)"
    assert source["row_factor_relation"] == "s^2=1-mass"
    assert source["base_normalization"] == "15/2 times the sectional Taylor coefficient before adding dual gain"
    pool = []
    for rows in source["polynomials"]:
        terms = {}
        for row in rows:
            assert len(row) == 3
            a, b = row[:2]
            assert type(a) is int and type(b) is int and min(a, b) >= 0
            assert (a, b) not in terms
            terms[(a, b)] = Q(row[2])
        pool.append(polynomial(terms))
    denominators, entries = {}, {}

    def entry(value):
        assert isinstance(value, list)
        if len(value) == 2 and all(type(i) is int for i in value):
            i, j = value
            assert 0 <= i < len(pool) and 0 <= j < len(pool)
            if j not in denominators:
                denominators[j] = circle_denominator(pool[j])
            if (i, j) not in entries:
                lead, a, b = denominators[j]
                entries[(i, j)] = RatCircle(pool[i]*fp(1/lead), a, b)
            return entries[(i, j)]
        return [entry(row) for row in value]

    data = {key: entry(source[key]) for key in ("A0", "A1", "d0", "d1", "baseR", "baseI")}
    for name in ("A0", "A1"):
        assert len(data[name]) == 20 and all(len(row) == 20 for row in data[name])
        assert all(isinstance(x, RatCircle) for row in data[name] for x in row)
    for name in ("d0", "d1"):
        assert len(data[name]) == 20 and all(isinstance(x, RatCircle) for x in data[name])
    for name in ("baseR", "baseI"):
        assert len(data[name]) == 3 and all(isinstance(x, RatCircle) for x in data[name])
    return data, {"polynomial_count": len(pool), "distinct_denominators_checked": len(denominators),
                  "distinct_rational_entries": len(entries)}


def dot(a, b):
    assert len(a) == len(b)
    value = RatCircle.scalar(0)
    for x, y in zip(a, b):
        if x and y:
            value = value+x*y
    return value


def matvec(matrix, vector, progress=lambda phase: None, phase="matrix-vector"):
    result = []
    for i, row in enumerate(matrix):
        result.append(dot(row, vector))
        progress(phase+"-row-"+str(i+1))
    return result


def quadratic(matrix, vector, progress=lambda phase: None, phase="quadratic"):
    assert len(matrix) == len(vector) and all(len(row) == len(vector) for row in matrix)
    # Generic complete quadratic form: all matrix rows, no special products.
    return dot(vector, matvec(matrix, vector, progress, phase))


def fraction_inverse(matrix, progress=lambda phase: None):
    n = len(matrix)
    assert n and all(len(row) == n for row in matrix)
    original = [[Q(x) for x in row] for row in matrix]
    work = [row+[Q(i == j) for j in range(n)] for i, row in enumerate(original)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if work[i][j]), None)
        assert pivot is not None, "constant bulk matrix must be invertible"
        work[j], work[pivot] = work[pivot], work[j]
        scale = work[j][j]
        work[j] = [x/scale for x in work[j]]
        for i in range(n):
            if i != j and work[i][j]:
                factor = work[i][j]
                work[i] = [x-factor*y for x, y in zip(work[i], work[j])]
        progress("Fraction-Gauss-Jordan-pivot-"+str(j+1))
    inverse = [row[n:] for row in work]
    for i in range(n):
        for j in range(n):
            assert work[i][j] == Q(i == j)
            assert sum((original[i][k]*inverse[k][j] for k in range(n)), ZERO) == Q(i == j)
            assert sum((inverse[i][k]*original[k][j] for k in range(n)), ZERO) == Q(i == j)
    return inverse


def antipodal_numerator(value, axis):
    exponent = (value.a, value.b)[axis]
    out = {}
    for mon, coefficient in coefficient_dict(value.n).items():
        power = mon[axis]
        assert power <= 2*exponent, "the angular function must extend continuously to the projective endpoint"
        transformed = list(mon)
        transformed[axis] = 2*exponent-power
        out[tuple(transformed)] = coefficient*((-1)**power)*Q(2)**(exponent-power)
    return polynomial(out)


def check_antipodal(value):
    for axis in (0, 1):
        assert antipodal_numerator(value, axis) == value.n, "all coefficients of the independent antipodal identity"


def verify_irrational_model(base_i):
    t = {0: ONE, 1: Q(16), 2: Q(-12), 3: Q(-32), 4: Q(4)}
    c = {0: ONE, 2: Q(-12), 4: Q(4)}
    s = {1: Q(4), 3: Q(-8)}
    functions = []
    for axis in (0, 1):
        def univariate(values):
            return polynomial({((k, 0) if axis == 0 else (0, k)): value for k, value in values.items()})
        numerator, cosine, sine = univariate(t), univariate(c), univariate(s)
        denominator = 4*circle_power(2 if axis == 0 else 0, 2 if axis == 1 else 0)
        assert numerator == cosine+4*sine
        assert cosine*cosine+2*sine*sine == denominator*denominator
        assert 9*denominator*denominator-numerator*numerator == 2*(2*cosine-sine)**2, "exact square proves |T| <= 3 at every real phase"
        functions.append(RatCircle(numerator*fp(Q(1, 4)), 2 if axis == 0 else 0, 2 if axis == 1 else 0))
    expected = Q(11, 160)*(3*functions[0]+7*functions[1])
    assert base_i[0] == Q(15, 2)*expected
    assert base_i[1] == -Q(15)*expected
    assert base_i[2] == 0
    return {"all_three_mass_coefficients_exact": True,
            "I": "11/160*(1-2m)*(3T(u)+7T(v))",
            "T": "(1+16u-12u^2-32u^3+4u^4)/(1+2u^2)^2",
            "T_square_identity_checked": 2, "uniform_absolute_I_bound": "33/16"}


def reconstruct_five(data, candidate, progress):
    reference = Q(candidate["reference_mass"])
    comparison = Q(candidate["sqrt2_rational_comparison"])
    margin = Q(candidate["rational_margin"])
    order = int(candidate["radial_bernstein_degree"])
    coefficients = [Q(x) for x in candidate["filter_coefficients"]]
    assert reference == Q(1, 3) and comparison == Q(7, 5) and margin == Q(1, 8)
    assert order == 4 and coefficients == [Q(17, 8), Q(-21, 16), Q(1, 4)]
    A0, A1, d0, d1 = (data[name] for name in ("A0", "A1", "d0", "d1"))
    for i in range(20):
        for j in range(20):
            assert A0[i][j] == A0[j][i] and A1[i][j] == A1[j][i]
            A1[i][j].as_fraction()
        progress("full-matrix-symmetry-and-constant-slope-row-"+str(i+1))
    Aref = [[A0[i][j]+reference*A1[i][j] for j in range(20)] for i in range(20)]
    dref = [d0[i]+reference*d1[i] for i in range(20)]
    bulk = [[Aref[i][j].as_fraction() for j in range(16)] for i in range(16)]
    inverse = fraction_inverse(bulk, progress)
    cross = [row[16:] for row in Aref[:16]]
    inv_cross = [[dot(inverse[i], [cross[j][k] for j in range(16)]) for k in range(4)] for i in range(16)]
    bulk_drive = matvec(inverse, dref[:16], progress, "independent-bulk-drive")
    small, reduced_drive = [], []
    for i in range(4):
        left = [cross[j][i] for j in range(16)]
        small.append([Aref[16+i][16+j]-dot(left, [inv_cross[k][j] for k in range(16)]) for j in range(4)])
        reduced_drive.append(dref[16+i]-dot(left, bulk_drive))
        progress("four-core-Schur-definition-row-"+str(i+1))
    assert all(small[i][j] == small[j][i] for i in range(4) for j in range(4))
    first_power = matvec(small, reduced_drive, progress, "S-times-b")
    second_power = matvec(small, first_power, progress, "S-squared-times-b")
    retained = [(34*reduced_drive[i]-21*first_power[i]+4*second_power[i])/16 for i in range(4)]
    remaining = [dref[i]-dot(cross[i], retained) for i in range(16)]
    x = matvec(inverse, remaining, progress, "frozen-witness-bulk")+retained
    assert len(x) == 20
    for i in range(16):
        assert dot(Aref[i][:16], x[:16])+dot(cross[i], x[16:]) == dref[i], "all sixteen frozen-reference bulk equations"
    # Rebuild the dual gain with the complete 20-dimensional quadratic forms.
    g0 = 2*dot(d0, x)-quadratic(A0, x, progress, "generic-A0-quadratic")
    g1 = 2*dot(d1, x)-quadratic(A1, x, progress, "generic-A1-quadratic")
    base = [data["baseR"][i]+comparison*data["baseI"][i] for i in range(3)]
    # x_actual=s*x, d_actual=s*(d0+m*d1), s^2=1-m.  No gain/slope shortcut.
    power = [base[0]-Q(15, 2)*margin+g0, base[1]+g1-g0, base[2]-g1]
    quadratic_bernstein = [power[0], power[0]+power[1]/2, sum(power, RatCircle.scalar(0))]
    elevated = [power[0]+Q(k, order)*power[1]+Q(k*(k-1), order*(order-1))*power[2]
                for k in range(order+1)]
    # A separate binomial formula checks the degree elevation itself.
    for k, value in enumerate(elevated):
        recovered = RatCircle.scalar(0)
        for j in range(3):
            if 0 <= k-j <= order-2:
                recovered += Q(math.comb(2, j)*math.comb(order-2, k-j), math.comb(order, k))*quadratic_bernstein[j]
        assert recovered == value
    return elevated, {"model_dimension": 20, "bulk_dimension": 16, "core_dimension": 4,
                      "constant_bulk_inverse_both_products_exact": True,
                      "generic_full_quadratic_forms_rebuilt": 2,
                      "frozen_bulk_equations_checked": 16,
                      "mass_power_coefficients_rebuilt": 3, "elevated_coefficients_rebuilt": len(elevated),
                      "filter": "(34I-21S+4S^2)b/16",
                      "gain_identity": "(1-m)*(2(d0+m*d1)^T x0-x0^T(A0+m*A1)x0)"}


@lru_cache(maxsize=512)
def direct_axis_transform(degree, lo, hi):
    """One closed formula combines direct affine expansion with binomial conversion."""
    assert lo < hi
    common = math.lcm(lo.denominator, hi.denominator)
    a = lo.numerator*(common//lo.denominator)
    b = hi.numerator*(common//hi.denominator)
    choose_scale = math.lcm(*(math.comb(degree, k) for k in range(degree+1)))
    rows = []
    for j in range(degree+1):
        rows.append([sum(math.comb(k, h)*a**(k-h)*(b-a)**h*common**(degree-k)
                         *math.comb(j, h)*(choose_scale//math.comb(degree, h))
                         for h in range(min(j, k)+1)) for k in range(degree+1)])
    return rows, common**degree*choose_scale


def direct_box_coefficients(power, box):
    degrees = tuple(max((mon[j] for mon in power), default=0) for j in (0, 1))
    du, dv = degrees
    original = [[0]*(dv+1) for _ in range(du+1)]
    for (a, b), c in power.items():
        assert type(c) is int
        original[a][b] = c
    uu, uscale = direct_axis_transform(du, *box[0])
    vv, vscale = direct_axis_transform(dv, *box[1])
    middle = [[sum(uu[i][a]*original[a][b] for a in range(du+1)) for b in range(dv+1)] for i in range(du+1)]
    answer = [[sum(middle[i][b]*vv[j][b] for b in range(dv+1)) for j in range(dv+1)] for i in range(du+1)]
    return answer, uscale*vscale


def verify_cover(certificate, powers, angular_box, progress):
    tree = certificate["tree"]
    assert isinstance(tree, list) and 0 < len(tree) <= 100000
    index, leaves, bounds, coefficient_count = 0, 0, 0, 0
    comparisons = []

    def visit(pending, box, depth):
        nonlocal index, leaves, bounds, coefficient_count
        assert depth <= 40 and index < len(tree)
        node = tree[index]
        index += 1
        closed, axis = node["closed"], node["axis"]
        assert isinstance(closed, list) and all(isinstance(x, str) for x in closed)
        assert len(set(closed)) == len(closed) and set(closed) <= pending
        assert type(axis) is int
        for name in closed:
            tensor, scale = direct_box_coefficients(powers[name], box)
            assert scale > 0 and all(c > 0 for row in tensor for c in row), "every claimed strict Bernstein coefficient must be positive after direct power reconstruction"
            du, dv = len(tensor)-1, len(tensor[0])-1
            minimum, i, j = min((tensor[i][j], i, j)
                                for i in range(du+1) for j in range(dv+1))
            count = (du+1)*(dv+1)
            coefficient_count += count
            comparisons.append({"node": index-1, "function": name,
                                "box": [[str(lo), str(hi)] for lo, hi in box],
                                "bidegree": [du, dv], "coefficient_count": count,
                                "minimum_Bernstein_coefficient": str(Q(minimum, scale)),
                                "minimum_location": [i, j]})
            bounds += 1
            progress("direct-Bernstein-bound-"+str(bounds))
        remaining = pending-set(closed)
        if axis == -1:
            assert not remaining, "every polynomial must be closed on every leaf"
            leaves += 1
            return
        assert axis in (0, 1) and remaining
        lo, hi = box[axis]
        middle = (lo+hi)/2
        for interval in ((lo, middle), (middle, hi)):
            child = list(box)
            child[axis] = interval
            visit(remaining, tuple(child), depth+1)

    visit(set(powers), angular_box, 0)
    assert index == len(tree), "the entire submitted tree is consumed, with both children of every split"
    assert index == 2*leaves-1, "a full single-root binary tree"
    return {"checked_boxes": index, "checked_leaf_boxes": leaves, "checked_polynomial_bounds": bounds,
            "checked_Bernstein_coefficients": coefficient_count,
            "comparisons": comparisons,
            "complete_tree_consumed": True, "all_strict_bounds_exact": True}


def self_tests():
    a = RatCircle(polynomial({(1, 0): Q(2), (0, 1): Q(-3), (0, 0): ONE}), 1, 2)
    b = RatCircle(polynomial({(2, 1): ONE, (0, 0): Q(-2)}), 2, 1)
    point = (Q(2, 5), Q(-3, 7))
    for value, expected in ((a+b, a.at(*point)+b.at(*point)),
                            (a*b, a.at(*point)*b.at(*point)),
                            (Q(3, 2)*a-b/3, Q(3, 2)*a.at(*point)-b.at(*point)/3)):
        assert value.at(*point) == expected
    denominator = 3*circle_power(2, 1)
    assert circle_denominator(denominator) == (Q(3), 2, 1)
    assert RatCircle(5*circle_power(2, 3), 2, 3).as_fraction() == 5
    fraction_inverse([[ZERO, Q(2), ONE], [ONE, ONE, ZERO], [Q(2), ZERO, Q(3)]])
    t = RatCircle(polynomial({(0, 0): Q(1, 4), (1, 0): Q(4), (2, 0): Q(-3),
                              (3, 0): Q(-8), (4, 0): ONE}), 2, 0)
    check_antipodal(t)
    noninvariant = RatCircle(polynomial({(1, 0): ONE}), 1, 0)
    assert antipodal_numerator(noninvariant, 0) != noninvariant.n
    assert antipodal_numerator(RatCircle(antipodal_numerator(noninvariant, 0), 1, 0), 0) == noninvariant.n
    power = {(0, 0): 3, (2, 1): 2, (1, 3): -4, (0, 2): 1}
    box = ((Q(-3, 4), Q(1, 8)), (Q(1, 4), Q(5, 8)))
    tensor, scale = direct_box_coefficients(power, box)
    du, dv = len(tensor)-1, len(tensor[0])-1
    for u, v in ((ZERO, ONE), (Q(1, 3), Q(2, 5))):
        reconstructed = sum((Q(tensor[i][j], scale)*math.comb(du, i)*u**i*(1-u)**(du-i)
                             *math.comb(dv, j)*v**j*(1-v)**(dv-j)
                             for i in range(du+1) for j in range(dv+1)), ZERO)
        physical = (box[0][0]+(box[0][1]-box[0][0])*u, box[1][0]+(box[1][1]-box[1][0])*v)
        expected = sum((c*physical[0]**i*physical[1]**j for (i, j), c in power.items()), ZERO)
        assert reconstructed == expected
    test_cover = {"tree": [{"axis": 0, "closed": []}, {"axis": -1, "closed": ["positive"]},
                           {"axis": -1, "closed": ["positive"]}]}
    test = verify_cover(test_cover, {"positive": {(0, 0): 2, (1, 0): 1}},
                        ((Q(-3, 4), Q(3, 4)), (Q(-3, 4), Q(3, 4))), lambda phase: None)
    assert test["checked_boxes"] == 3 and test["checked_polynomial_bounds"] == 2
    rejected = 0
    for bad_tree, bad_power in (
            ({"tree": test_cover["tree"][:-1]}, {"positive": {(0, 0): 2}}),
            ({"tree": [{"axis": -1, "closed": []}]}, {"positive": {(0, 0): 2}}),
            ({"tree": [{"axis": -1, "closed": ["positive"]}]}, {"positive": {(0, 0): -1}})):
        try:
            verify_cover(bad_tree, bad_power, box, lambda phase: None)
        except AssertionError:
            rejected += 1
    assert rejected == 3, "truncated covers, uncovered functions, and negative bounds must fail"
    return {"RatCircle_arithmetic_evaluations": 3, "circle_denominator_recognition": 1,
            "constant_recognition": 1, "Fraction_inverse_with_row_pivot": 1,
            "antipodal_positive_negative_and_involution_checks": 3,
            "direct_affine_Bernstein_evaluations": 2, "complete_tree_check": 1,
            "invalid_certificates_rejected": rejected}


def primitive_numerator(value):
    """Derive the integer numerator; no saved polynomial is an input."""
    assert value.n, "an identically zero coefficient cannot prove strict positivity"
    denominator = circle_power(value.a, value.b)
    common = value.n.gcd(denominator)
    numerator, denominator = value.n/common, denominator/common
    lead, a, b = circle_denominator(denominator)
    assert lead > 0
    coefficients = coefficient_dict(numerator)
    scale = math.lcm(*(c.denominator for c in coefficients.values()))
    integers = {mon: c.numerator*(scale//c.denominator) for mon, c in coefficients.items()}
    content = math.gcd(*map(abs, integers.values()))
    assert content > 0
    integers = {mon: c//content for mon, c in integers.items()}
    rebuilt = RatCircle(polynomial(integers)*fp(Q(content, scale)/lead), a, b)
    assert value.cross_equal(rebuilt), "the derived primitive polynomial reproduces the full witness"
    assert math.gcd(*map(abs, integers.values())) == 1
    check_antipodal(value)
    check_antipodal(rebuilt)
    return integers, (a, b)


def check_b_bound():
    """Check the two elementary source norms in the direct B hand proof."""
    norm_residual = Q(49, 384)+Q(27, 128)+Q(25, 24)+Q(27, 64)+Q(27, 16000)
    norm_vertical = Q(625, 384)+Q(27, 128)+Q(2, 3)+Q(27, 64)
    assert norm_residual == Q(173, 96)+Q(27, 16000)
    assert norm_vertical == Q(281, 96)
    cost = 2*norm_residual+Q(2, 3)*norm_vertical
    assert cost == Q(50, 9)+Q(27, 8000) and cost < Q(45, 8)
    lower = Q(1, 8)-cost/54
    assert lower > Q(1, 48)
    return {"dual_cost": str(cost), "second_derivative_lower_bound": str(lower),
            "exceeds_one_forty_eighth": True,
            "geometric_source_identity": "the direct B argument in the accompanying proof"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="optional JSON result file")
    parser.add_argument("--status", type=Path, help="optional live JSON progress file")
    args = parser.parse_args()
    if not __debug__:
        parser.error("exact verification requires assertions; do not run with -O")
    began, last_save = time.monotonic(), [0.0]
    for path in (args.output, args.status):
        if path:
            path.parent.mkdir(parents=True, exist_ok=True)
    result = {"status": "running", "canary": "GM-portable-geometry-to-sign-certificate-v1",
              "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, "gm_geometry.py")},
              "python": platform.python_version(), "flint": flint.__version__,
              "completed_steps": 0, "angular_checks": [], "external_data_files": 0,
              "scope": "rebuild the complete source model from metric formulas, then independently verify the final exact algebraic inequalities"}

    def save(phase, terminal=False):
        event = {"status": result["status"], "phase": phase, "completed": result["completed_steps"],
                 "total": 4, "heartbeat": datetime.now(timezone.utc).isoformat(),
                 "elapsed_seconds": time.monotonic()-began, "terminal": terminal,
                 "expected_finish": "unknown"}
        for path, value in ((args.output, result), (args.status, event)):
            if path:
                temporary = path.with_suffix(path.suffix+".tmp")
                temporary.write_text(json.dumps(value, indent=2), encoding="utf-8")
                temporary.replace(path)
        print(f"[{result['completed_steps']}/4] {phase} ({time.monotonic()-began:.1f}s)", flush=True)
        last_save[0] = time.monotonic()

    def progress(phase):
        if time.monotonic()-last_save[0] > 10:
            save(phase)

    save("Starting exact verification")
    try:
        from gm_geometry import build_model
        import numpy
        import sympy
        result["numpy"] = numpy.__version__
        result["sympy"] = sympy.__version__
        result["self_tests"] = self_tests()
        save("Rebuilding the source geometry and all mass identities")
        source = build_model(progress)
        result["geometry_checks"] = source.get("exact_checks", {})
        result["geometry_source_provenance"] = source.get("source_provenance", [])
        # The geometry constructor returns coefficient data in memory. There
        # are no saved models, angular numerators, or producer success flags.
        data, statistics = decode_model(source)
        result["model_statistics"] = statistics
        result["irrational_model"] = verify_irrational_model(data["baseI"])
        result["completed_steps"] = 1
        save("Geometry and irrational-part identity verified")
        exact, reconstruction = reconstruct_five(data, WITNESS, progress)
        result["model_reconstruction"] = reconstruction
        result["completed_steps"] = 2
        save("Five angular functions independently reconstructed")
        bx = tuple(tuple(Q(x) for x in pair) for pair in WITNESS["angular_box"])
        assert bx == ((Q(-3, 4), Q(3, 4)), (Q(-3, 4), Q(3, 4)))
        assert all(lo == -hi and hi*hi > Q(1, 2) for lo, hi in bx)
        assert len(exact) == 5
        powers = {}
        for index, value in enumerate(exact):
            name = "mass-coefficient-"+str(index)
            coefficients, exponents = primitive_numerator(value)
            degrees = [max(mon[j] for mon in coefficients) for j in (0, 1)]
            assert degrees == ([32, 36] if index < 4 else [8, 8]), "the primitive bidegrees must agree with the paper"
            assert (degrees[0]+1)*(degrees[1]+1) < 25000
            powers[name] = coefficients
            result["angular_checks"].append({"name": name, "degrees": degrees, "terms": len(coefficients),
                                              "circle_exponents": list(exponents),
                                              "primitive_prefactor_identity_exact": True,
                                              "both_antipodal_identities_exact": True,
                                              "strictly_positive_circle_denominator": True})
            progress("Reconstructing integer numerator "+str(index+1)+"/5")
        result["projective_cover"] = {"box": WITNESS["angular_box"], "squared_halfwidth": "9/16",
                                       "antipodal_threshold": "1/2", "proper_endpoint_extensions_checked": True,
                                       "reason": "outside the square, each needed involution u -> -1/(2u) sends that coordinate inside it; infinity maps to zero"}
        certificate = {"tree": [{"axis": axis, "closed": ["mass-coefficient-"+str(i) for i in indices]}
                                for axis, indices in PARTITION]}
        result["box_verification"] = verify_cover(certificate, powers, bx, progress)
        cover = result["box_verification"]
        assert (cover["checked_boxes"], cover["checked_leaf_boxes"], cover["checked_polynomial_bounds"]) == (29, 15, 38), "the verified subdivision must agree with the paper"
        result["completed_steps"] = 3
        save("All 38 exact sign bounds and the complete angular cover verified")
        lower, upper = Q(7, 5), Q(10, 7)
        absolute_i = Q(11, 160)*(3+7)*3
        assert absolute_i == Q(33, 16)
        assert lower*lower < 2 < upper*upper
        derived = Q(WITNESS["rational_margin"])-(upper-lower)*absolute_i
        assert derived == Q(37, 560) and derived > Q(1, 16)
        assert 2*derived > Q(1, 8)
        result["B_bound"] = check_b_bound()
        result.update({"status": "complete", "completed_steps": 4,
                       "uniform_rational_witness_margin_verified": "1/8",
                       "sqrt2_interval": [str(lower), str(upper)],
                       "derived_A_Taylor_coefficient_lower_bound": str(derived),
                       "A_Taylor_coefficient_exceeds_one_sixteenth": True,
                       "portable_certificate_verified": True,
                       "elapsed_seconds": time.monotonic()-began,
                       "decision": "The source formulas generate the global mass model; independent reconstruction and complete exact sign checks prove a_A > 37/560. The direct B estimate gives reduced second derivative > 1/48."})
        save("PASS: exact algebraic verification complete", True)
        print("A: Taylor coefficient > 37/560; reduced second derivative > 1/8.")
        print("B: reduced second derivative > 1/48.")
    except Exception as error:
        result.update({"status": "failed", "error": repr(error), "elapsed_seconds": time.monotonic()-began})
        save("FAIL: "+str(error), True)
        raise


if __name__ == "__main__":
    main()
