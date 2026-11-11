# Code and data for "A lower bound of order A^(4/3) on the number of mass points of the capacity-achieving input of the amplitude-constrained Gaussian channel"

This repository holds every computation that the paper mentions. The paper itself is distributed separately.

No step of the proof of Theorem 1 uses a computer. The scripts in `proofs/` check numbers that the text derives by hand, and they give the sharper constant of Appendix C.

## Contents
- `proofs/check_hand_J2.py`: checks the numbers of the hand proof of Lemma 2(c), the depth 2. It covers the closed form of the bound J, the three hand inequalities, the constants of Theorem 1 and the chord bound on random Gaussian mixtures. Its output is `check_hand_J2_out.txt`.
- `proofs/sharp_depth_ball_arithmetic.py`: the computation of Appendix C in ball arithmetic. It shows that J(1.1648) < 0 and J(1.1647) > 0, and gives the constant 0.1660. Its output is `sharp_depth_ball_arithmetic_out.txt`.
- `computations/data/`: the 2668 computed optimal inputs of Appendix D (`stored-optima/`, one file each) and the index `count_law_data.json`, which lists them with their amplitudes and counts.
- `computations/scripts/`: the scripts of Appendix D and of the figures.
  - `summary_states.py` writes `numbers_states.json`: the number of stored optima, their largest Kuhn-Tucker residual, and the ratio of the computed counts to the proved bound.
  - `numbers_A2.py` writes `numbers_A2.json`: the parts of J(1.1648) and the differences between the computed counts and 0.4002 (A + beta/2)^(4/3).
  - `make_fig_face.py`, `make_fig_step2.py`, `make_fig_flat.py` and `make_fig_count.py` draw Figs. 1 to 4 of the paper, each with the numbers of its caption.

## How to run
Python 3 with the packages of `requirements.txt` (python-flint for ball arithmetic, mpmath, numpy, matplotlib).

    python proofs/check_hand_J2.py
    python proofs/sharp_depth_ball_arithmetic.py
    python computations/scripts/summary_states.py

The last line of `check_hand_J2.py` reads `CHECK hand_J2: ALL FOUR STATEMENTS HOLD`, and that of `sharp_depth_ball_arithmetic.py` reads `CHECK sharp_depth: ALL THREE STATEMENTS HOLD`. The scripts in `computations/scripts/` write their outputs next to themselves.

## License
See `LICENSE`.
