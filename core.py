"""Shared Slipway algebra, search solvers, and certificate checks."""

from time import perf_counter

from sage.all import GF, PolynomialRing, is_prime, matrix, vector

from reference import Poseidon


def require(condition, message):
    """Keep certificate checks enabled even when Python runs with -O."""
    if not condition:
        raise ValueError(message)


def validate_parameters(p, t, rf, rp):
    require(p > 2 * t and is_prime(p) and p % 3 == 2,
            "Require a prime p > 2*width with p = 2 mod 3 (bijective cubing).")
    require(rf >= 4 and rf % 2 == 0 and rp >= 0,
            "Require even rf >= 4 and rp >= 0.")
    require(t >= 3 * (rf // 2), "Prefix construction requires width >= 3*rf/2.")


class Slipway:
    """Cubic CICO-2 with a homogeneous prefix and a maximal finite trail."""

    def __init__(self, data):
        self.p, self.t = data["prime"], data["width"]
        self.rf, self.rp = data["full_rounds"], data["partial_rounds"]
        validate_parameters(self.p, self.t, self.rf, self.rp)
        self.field = GF(self.p)
        self.matrix = matrix(self.field, data["matrix"])
        require(self.matrix.dimensions() == (self.t, self.t), "Wrong matrix shape.")
        require(len(data["round_constants"]) == self.rf + self.rp
                and all(len(row) == self.t for row in data["round_constants"]),
                "Wrong round-constant shape.")
        self.constants = [vector(self.field, row) for row in data["round_constants"]]
        self.d_c, self.d_x, self.d_w, self.u1, self.u2 = (
            vector(self.field, data[key]) for key in ("d_c", "d_x", "d_w", "u1", "u2")
        )
        self.ratios = [[self.field(v) for v in pair] for pair in data["ratios"]]
        require(len(self.ratios) == self.rf // 2 - 1, "Wrong prefix length.")
        self.control = data.get("control")
        # New searches supply upstream-generated constants; the paper fixture
        # supplies its reported constants. Both use the supplied matrix.
        self.reference = Poseidon(
            prime=self.p, alpha=3, t=self.t, r_f=self.rf, r_p=self.rp,
            mds=[[int(v) for v in row] for row in self.matrix.rows()],
            round_constants=[int(c) for row in self.constants for c in row],
        )

    def direction(self, control):
        """Coefficient of X after the initial full rounds."""
        g, h = 1, control
        for lambda_g, lambda_h in self.ratios:
            g, h = (g + lambda_g * h)**3, (g + lambda_h * h)**3
        return g * self.u1 + (h - g) * self.u2

    def round(self, state, index):
        """Add constants, cube all/one coordinate, then apply the matrix."""
        state = state + vector(state.base_ring(), self.constants[index])
        active = range(self.t) if index < self.rf // 2 or index >= self.rf // 2 + self.rp else (0,)
        for i in active:
            state[i] = state[i]**3
        return self.matrix * state

    def recover_input(self, root, control):
        """Invert X=x^(3^(RF/2-1)) and the first S-box, both bijective."""
        x = self.field(root)**pow(3**(self.rf // 2 - 1), -1, self.p - 1)
        q0 = self.d_c + x * (self.d_x + self.field(control) * self.d_w)
        inverse_cube = pow(3, -1, self.p - 1)
        return [int(v**inverse_cube - c) for v, c in zip(q0, self.constants[0])]

    def certify(self):
        """Exact prefix identity, trail dimension, and power irreducibility.

        These checks do not certify MDS: exhaustive minor enumeration is omitted.
        """
        ring = PolynomialRing(self.field, names=("x", "w"))
        x, w = ring.gens()
        q0 = vector(ring, self.d_c) + x * (self.d_x + w * self.d_w)
        require(all(q0[i] == self.constants[0][i]**3 for i in (0, 1)),
                "The family does not impose both input constraints.")
        state = self.matrix * q0  # q0 is already after the first S-box.
        for index in range(1, self.rf // 2):
            state = self.round(state, index)
        require(state == x**(3**(self.rf // 2 - 1)) * self.direction(w),
                "Prefix identity failed.")

        row = vector(self.field, [1] + [0] * (self.t - 1))
        rows = []
        for _ in range(self.t):
            rows.append(row)
            row = row * self.matrix
        observability = matrix(self.field, rows)
        require(observability.rank() == self.t, "Nonmaximal partial-round rank.")
        require(matrix(self.field, [self.u1, self.u2]).rank() == 2,
                "Dependent post-prefix directions.")
        require(observability[:self.t - 1] * self.u1 == 0
                and observability[:self.t - 2] * self.u2 == 0
                and self.matrix * self.u1 == self.u2, "Trail certificate failed.")
        power = self.matrix
        for exponent in range(1, 4 * self.t + 1):
            require(power.charpoly().is_irreducible(),
                    f"Reducible characteristic polynomial for M^{exponent}.")
            power *= self.matrix

    def output_polynomials(self, control):
        """Expand the actual tail, recording every partial S-box input degree."""
        ring = PolynomialRing(self.field, "X")
        state = vector(ring, self.direction(self.field(control))) * ring.gen()
        active_degrees = []
        for index in range(self.rf // 2, self.rf + self.rp):
            if index < self.rf // 2 + self.rp:
                active_degrees.append(max(0, state[0].degree()))
            state = self.round(state, index)
        return state[0], state[1], active_degrees

    def reference_permutation(self, inputs):
        """Evaluate every round with the unmodified poseidon-tools dependency."""
        return self.reference.permutation([int(v) for v in inputs])

    def recover_control(self, a, b):
        """Invert the prefix on a*u1+b*u2, returning (X,w) when w is finite."""
        g, h = self.field(a), self.field(a + b)
        inverse_cube = pow(3, -1, self.p - 1)
        for lg, lh in reversed(self.ratios):
            g, h = g**inverse_cube, h**inverse_cube
            g, h = (lh*g - lg*h)/(lh-lg), (h-g)/(lh-lg)
        if not g:
            return (0, 0) if not h else None  # The omitted projective direction.
        return int(g**(3**len(self.ratios))), int(h/g)

    def solve_joint(self, max_degree=81):
        """Solve both post-prefix coefficients together for a short active tail.

        A resultant eliminates A from the two output equations in A,B. Work
        in GF(p)[B][A] to use PARI, including at KoalaBear's characteristic.
        Every recovered control is checked by the fixed-control solver below.
        """
        bound = 3**(self.rf // 2 + max(0, self.rp - (self.t - 2)))
        require(bound <= max_degree,
                f"Joint degree bound {bound} exceeds --max-degree {max_degree}.")
        start = perf_counter()
        rb = PolynomialRing(self.field, "B")
        ra = PolynomialRing(rb, "A")
        state = ra.gen()*vector(ra, self.u1) + rb.gen()*vector(ra, self.u2)
        for index in range(self.rf // 2, self.rf + self.rp):
            state = self.round(state, index)
        f, g = state[:2]
        resultant = f.resultant(g)
        if not resultant:
            return None  # Degenerate elimination: retry with another matrix.
        univariate = PolynomialRing(self.field, "T")
        for b in sorted(resultant.roots(multiplicities=False)):
            common = univariate([c(b) for c in f]).gcd(univariate([c(b) for c in g]))
            roots = common.roots(multiplicities=False) if common else [self.field.zero()]
            for a in sorted(roots):
                recovered = self.recover_control(a, b)
                if recovered is None:
                    continue
                X, w = recovered
                require(X*self.direction(self.field(w)) == a*self.u1 + b*self.u2,
                        "Inverting the prefix failed.")
                result = self.solve(w)
                require(any(s["root"] == X for s in result["solutions"]),
                        "Joint root failed full-permutation verification.")
                result["resultant_degree"] = int(resultant.degree())
                result["joint_seconds"] = perf_counter() - start
                return result
        return None

    def solve(self, control):
        """Solve one fixed-control fiber; no enumeration of control values."""
        require(0 <= control < self.p, "Control must be a canonical field element.")
        start = perf_counter()
        f0, f1, active_degrees = self.output_polynomials(control)
        inactive = next((i for i, degree in enumerate(active_degrees) if degree), self.rp)
        expected = 3**(self.rf // 2 + self.rp - inactive)
        require(max(f0.degree(), f1.degree()) <= expected, "Degree bound exceeded.")
        polynomial_seconds = perf_counter() - start

        # Cross-check the polynomial representation against the unskipped evaluator.
        for root in (0, 1, self.p - 1):
            output = self.reference_permutation(self.recover_input(root, control))
            require([int(f0(root)), int(f1(root))] == output[:2],
                    "Polynomial and direct evaluations disagree.")

        start = perf_counter()
        X = f0.parent().gen()
        # Modular exponentiation never materializes the degree-p polynomial X^p.
        # Either output can be the zero polynomial for small parameters.
        polynomial = f0 if f0 else f1
        require(bool(polynomial), "Both output polynomials vanish identically.")
        root_factor = polynomial.gcd(pow(X, self.p, polynomial) - X)
        roots = sorted(int(root) for root in root_factor.roots(multiplicities=False))
        root_seconds = perf_counter() - start
        candidates, solutions = [], []
        for root in roots:
            inputs = self.recover_input(root, control)
            outputs = self.reference_permutation(inputs)
            require(inputs[:2] == [0, 0] and outputs[0] == 0, "Invalid root candidate.")
            require(outputs[1] == int(f1(root)), "Second-output evaluation disagrees.")
            candidates.append({"root": root, "outputs": outputs[:2]})
            if outputs[1] == 0:
                solutions.append({"root": root, "input": inputs, "output": outputs})
        return {
            "control": control,
            "inactive_partial_rounds": inactive,
            "output_degrees": [int(f.degree()) if f else None for f in (f0, f1)],
            "candidates": candidates,
            "solutions": solutions,
            "polynomial_seconds": polynomial_seconds,
            "root_seconds": root_seconds,
        }
