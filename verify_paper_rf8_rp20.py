"""Verify the supplied paper witness at KoalaBear, width 16, RF=8, RP=20.

Run with: sage -python verify_paper_rf8_rp20.py
The matrix, constants, control, root, input, and output are supplied data.
This verifier performs no matrix search, control search, or root extraction.
"""

import argparse
import json
from pathlib import Path
from time import perf_counter

from core import Slipway, require


INSTANCE = Path(__file__).with_name("paper_rf8_rp20.json")


def verify_paper(data):
    """Check the supplied paper witness and its algebraic certificates without solving."""
    start = perf_counter()
    instance = Slipway(data)
    require((instance.p, instance.t, instance.rf, instance.rp) == (2130706433, 16, 8, 20),
            "Expected the paper's KoalaBear RF=8, RP=20, width-16 instance.")
    witness = data["witness"]
    control, root = data["control"], witness["root"]
    inputs, outputs = witness["input"], witness["output"]
    require(0 <= control < instance.p and 0 <= root < instance.p,
            "Control and root must be canonical field elements.")
    require(all(len(values) == instance.t and all(0 <= v < instance.p for v in values)
                for values in (inputs, outputs)), "Invalid witness vector.")
    require(inputs[:2] == outputs[:2] == [0, 0], "Witness violates CICO-2 constraints.")
    require(instance.recover_input(root, control) == inputs,
            "Supplied root and control do not match the supplied input.")
    require(instance.reference_permutation(inputs) == outputs,
            "Supplied output does not match the reference permutation.")

    instance.certify()
    f0, f1, degrees = instance.output_polynomials(control)
    require(degrees[:14] == [0]*14 and degrees[14] == 1,
            "Expected fourteen inactive partial S-boxes.")
    require(f0.degree() == f1.degree() == 59049, "Expected output degree 59049.")
    require(f0(root) == f1(root) == 0, "Supplied root fails the output equations.")
    return {"control": control, "root": root, "seconds": perf_counter() - start}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    try:
        result = verify_paper(json.loads(INSTANCE.read_text()))
    except (OSError, KeyError, TypeError, ValueError) as error:
        parser.exit(1, f"Verification failed: {error}\n")
    print("VERIFIED: supplied paper witness; KoalaBear, t=16, RF=8, RP=20.")
    print(f"Supplied control: {result['control']}; supplied root: {result['root']}.")
    print("PASS: full input/output equality using poseidon-tools; CICO-2 constraints.")
    print("PASS: symbolic prefix; dim K14=2; powers 1..64 irreducible; output degrees 59049.")
    print("No matrix/control/root search. The 601080389 MDS minors are not rechecked.")
    print(f"Verification time after Sage import: {result['seconds']:.3f}s")


if __name__ == "__main__":
    main()
