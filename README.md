# Code and data for "A lower bound of order A^(4/3) on the number of mass points of the capacity-achieving input of the amplitude-constrained Gaussian channel"

This repository holds every computation that the paper mentions. The paper itself is distributed separately.

## Contents
- `proofs/check_hand_J2.py`: checks the numbers of the hand proof of Lemma 2(c), the depth 2. It covers the closed form of the bound J, the three hand inequalities, the constants of Theorem 1 and the chord bound on random Gaussian mixtures. Its output is `check_hand_J2_out.txt`.
- `proofs/q1796_lower_bound_constants_arb.py`: the computation of Appendix C in ball arithmetic. It shows that J(1.1648) < 0 and J(1.1647) > 0, and gives the constant 0.1660. Its output is `q1796_out.txt`.
- `computations/`: the computed optima of Appendix D (`data/`) and the scripts that draw the figures and compute the numbers of Appendix D (`scripts/`).

No step of the proof of Theorem 1 uses a computer. The scripts in `proofs/` check numbers that the text derives by hand, and they give the sharper constant of Appendix C.

## How to run
Python 3 with the packages of `requirements.txt` (python-flint for ball arithmetic, mpmath, numpy, matplotlib).

    py proofs/check_hand_J2.py
    py proofs/q1796_lower_bound_constants_arb.py

Each script prints its result. The last line of `check_hand_J2.py` reads `CHECK hand_J2: ALL FOUR STATEMENTS HOLD`.

## License
See `LICENSE`.
