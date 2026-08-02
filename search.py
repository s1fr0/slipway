"""Search a fresh Slipway matrix and CICO-2 control with configurable rounds.

Run with: sage -python search.py --rf 4 --rp 5
The default is a reduced-round KoalaBear instance; see README.md for scope.
"""

import argparse
from itertools import combinations
import json
from pathlib import Path
from random import Random
from time import perf_counter

from sage.all import GF, binomial, matrix, vector

from core import Slipway, require, validate_parameters
from reference import Poseidon, REVISION


def is_mds(M):
    """Check every nonempty square minor (923 at width six)."""
    for size in range(1, M.nrows() + 1):
        subsets = list(combinations(range(M.nrows()), size))
        for rows in subsets:
            for columns in subsets:
                if not M.matrix_from_rows_and_columns(rows, columns).det():
                    return False
    return True


def construct(rng, p, t, rf, rp, constants):
    """Absorb rf/2 initial full rounds, then complete a maximal finite trail.

    Each intermediate S-box separates G,H from two constant coordinates.
    The last stage sets u1=z_H and prescribes M(z_G+z_H)=u1, Mu1=u2.
    The prefix uses 3*rf/2-1 matrix images; the remaining images complete
    the chain u1,...,u_(t-1), all zero in coordinate 0.
    """
    F, half = GF(p), rf // 2

    def random_vector(nonzero=False):
        return vector(F, [rng.randrange(int(nonzero), p) for _ in range(t)])

    dc, dx, dw = random_vector(), random_vector(), random_vector()
    for i in (0, 1):
        dc[i], dx[i], dw[i] = F(constants[0][i])**3, 0, 0
    stages, ratios = [], []
    for j in range(half - 1):
        # Disjoint constant supports avoid dependent domains at longer prefixes.
        holes = () if j == half - 2 else (t - 2 - 2*j, t - 1 - 2*j)
        support = [i for i in range(t) if i not in holes]
        G, H = support[:len(support)//2], support[len(support)//2:]
        a = vector(F, [rng.randrange(1, p) if i in support else 0 for i in range(t)])
        lg, lh = [F(v) for v in rng.sample(range(1, p), 2)]
        b = vector(F, [v*(lg if i in G else lh) for i, v in enumerate(a)])
        e = vector(F, [rng.randrange(1, p) if i in holes else 0 for i in range(t)])
        zg = vector(F, [v**3 if i in G else 0 for i, v in enumerate(a)])
        zh = vector(F, [v**3 if i in H else 0 for i, v in enumerate(a)])
        stages.append((a, b, e, zg, zh))
        ratios.append([int(lg), int(lh)])

    domains = [dx, dw, dc]
    a, b, e, _, _ = stages[0]
    images = [a, b, e - vector(F, constants[1])]
    for j in range(len(stages) - 1):
        _, _, e, zg, zh = stages[j]
        a_next, b_next, e_next, _, _ = stages[j + 1]
        domains.extend([vector(F, [v**3 for v in e]), zg, zh])
        images.extend([e_next - vector(F, constants[j + 2]), a_next, b_next])
    _, _, _, zg, zh = stages[-1]
    chain = [zh]
    for _ in range(t - 3*half + 1):
        u = random_vector(nonzero=True)
        u[0] = 0
        chain.append(u)
    domains.extend([zg + zh, *chain])
    images.extend(chain)
    basis = matrix(F, domains).transpose()
    if not basis.is_invertible():
        return None
    inverse = basis.inverse()
    M0 = matrix(F, images + [vector(F, t)]).transpose() * inverse
    row, rows = vector(F, [1] + [0]*(t - 1)), []
    for _ in range(t - 1 - len(chain)):
        rows.append(row)
        row *= M0
    kernel = matrix(F, rows).right_kernel().basis()
    last = sum((F(rng.randrange(p))*v for v in kernel), vector(F, t))
    M = matrix(F, images + [last]).transpose() * inverse
    if not M.charpoly().is_irreducible() or not is_mds(M):
        return None
    data = {
        "prime": p, "width": t, "full_rounds": rf, "partial_rounds": rp,
        "matrix": [[int(v) for v in row] for row in M.rows()],
        "round_constants": constants, "ratios": ratios,
        "constants_source": f"khovratovich/poseidon-tools@{REVISION}",
    }
    for name, value in (("d_c", dc), ("d_x", dx), ("d_w", dw),
                        ("u1", chain[0]), ("u2", chain[1])):
        data[name] = [int(v) for v in value]
    instance = Slipway(data)
    try:
        instance.certify()
    except ValueError:
        return None
    return data, instance


def search(seed=0, max_attempts=10000, *, p=2130706433, t=None, rf=4, rp=5,
           max_degree=81):
    """Fix upstream constants, search matrices, and solve the two output equations."""
    start = perf_counter()
    t = 3*(rf // 2) if t is None else t
    validate_parameters(p, t, rf, rp)
    require(max_attempts > 0, "--max-attempts must be positive.")
    bound = 3**(rf // 2 + max(0, rp - (t - 2)))
    require(bound <= max_degree,
            f"Joint degree bound {bound} exceeds --max-degree {max_degree}.")
    # Generated once, before any matrix is sampled; no local constant generator.
    constants = Poseidon(p, 3, t, rf, rp).round_constants
    rng = Random(seed)
    matrices = 0
    for attempt in range(1, max_attempts + 1):
        candidate = construct(rng, p, t, rf, rp, constants)
        if candidate is None:
            continue
        data, instance = candidate
        matrices += 1
        result = instance.solve_joint(max_degree)
        if result is not None:
            data["control"] = result["control"]
            return {"instance": data, "result": result, "seed": seed,
                    "matrix_attempts": attempt, "certified_matrices": matrices,
                    "seconds": perf_counter() - start}
    raise RuntimeError(f"No solution within {max_attempts} matrix attempts "
                       f"({matrices} certified matrices).")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime", type=int, default=2130706433, help="Prime field; default KoalaBear.")
    parser.add_argument("--width", type=int, help="State width; default 3*rf/2.")
    parser.add_argument("--rf", type=int, default=4, help="Even number of full rounds, at least 4.")
    parser.add_argument("--rp", type=int, default=5, help="Number of partial rounds.")
    parser.add_argument("--seed", type=int, default=0, help="Reproducible matrix search seed.")
    parser.add_argument("--max-attempts", type=int, default=10000, help="Bound matrix attempts.")
    parser.add_argument("--max-degree", type=int, default=81, help="Bound joint-solver polynomial degree.")
    parser.add_argument("--output", type=Path, help="Save the instance, solution, and search counts.")
    args = parser.parse_args()
    try:
        report = search(args.seed, args.max_attempts, p=args.prime, t=args.width,
                        rf=args.rf, rp=args.rp, max_degree=args.max_degree)
    except (RuntimeError, ValueError) as error:
        parser.exit(1, f"{error}\n")
    data, result = report["instance"], report["result"]
    t = data["width"]
    print(f"Fresh GF({data['prime']}), t={t}, alpha=3, RF={args.rf}, RP={args.rp} instance.")
    print(f"PASS: all {binomial(2*t, t)-1} minors; symbolic prefix; "
          f"dim K{t-2}=2; powers 1..{4*t} irreducible.")
    print(f"Matrix attempts: {report['matrix_attempts']}; certified matrices: "
          f"{report['certified_matrices']}; resultant degree: {result['resultant_degree']}")
    print(f"Control: {result['control']}; inactive partial rounds: "
          f"{result['inactive_partial_rounds']}; degrees: {result['output_degrees']}")
    for solution in result["solutions"]:
        print(f"CICO-2 solution X={solution['root']}")
        print(f"Input:  {solution['input']}")
        print(f"Output: {solution['output']}")
    print(f"Search time after Sage import: {report['seconds']:.3f}s")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
