# Two details that decide whether a port of airGR is identical

Translating the published GR equations gets you a model that behaves correctly.
It does not get you airGR's numbers. Two implementation details, invisible in the
papers, stand between the two — and both were found the same way: by staring at a
deviation that was far too structured to be noise.

If you are porting airGR to another language, this page is the part worth reading.

---

## 1. Fortran evaluates `0.9` in single precision

airGR's kernel splits effective rainfall between its two routing branches using a
constant declared like this:

```fortran
doubleprecision, parameter :: B=0.9
```

The literal `0.9` carries no precision suffix, so the compiler evaluates it as a
default-kind `REAL(4)` and *then* promotes the result to double. `B` does not hold
0.9. It holds:

```
0.89999997615814208984375
```

A port that writes `B = 0.9` in a language with double-precision literals is not
reproducing airGR — it is reproducing what the airGR authors meant, which is a
different thing. The discharge series then drifts by a relative deviation of

```
2.384185791015625e-07  =  2**-22
```

That exponent is the tell. `2⁻²²` is the epsilon of single precision; seeing it
appear in a double-precision comparison is a signature, not a coincidence, and it
is what pointed at the cause here.

### Where it applies

| Constant | Model | Fortran source | Actual double value |
|---|---|---|---|
| `B = 0.9` | GR4J, GR5J, GR6J, GR4H, GR5H | routing split | 0.899999976158142 |
| `C = 0.4` | GR6J | exponential store share | 0.400000005960464 |
| `1./3.` | GR2M | percolation exponent | 0.333333343267441 |
| `0.7`, `0.3` | GR1A | rainfall weighting | 0.699999988079071, 0.300000011920929 |
| `MinSpeed = 0.1` | CemaNeige | melt rate floor | 0.100000001490116 |
| `0.9 * MeanAnSolidPrecip` | CemaNeige | melt threshold | — |
| `-999.999` | CemaNeige | sentinel for `Glocalmax` | −999.9990234375 |

That last row is how the snow module gave itself away. airGR's own output files
contain `Glocalmax = -999.9990234375` — the single-precision image of `-999.999`,
printed in a double-precision column.

Constants that *are* exactly representable in `REAL(4)` — `2.5`, `13.`, `25.62890625`
(which is `(9/4)⁴`), `759.69140625` (`(21/4)⁴`) — pass through unharmed. Only the
ones that are not, matter.

### How GRsuite handles it

```python
_F32_09 = float(np.float32(0.9))   # 0.8999999761581421
```

Each affected constant is promoted once, at import, with a comment saying why.
See [`src/grsuite/_kernels.py`](../src/grsuite/_kernels.py).

### If you are porting airGR yourself

Grep the Fortran for real literals without a `_8`, `d0` or `_dp` suffix that end up
in double-precision arithmetic. Every one of them is a landmine. The check that
catches this class of bug in one pass: run your port against airGR and look at the
*exponent* of the relative deviation. If it is near `10⁻⁷`, you have a single-precision
literal somewhere. If it is near `10⁻¹⁶`, you are done.

---

## 2. The automatic warm-up loses a day to leap years

When the user does not specify a warm-up period, airGR builds one by stepping back
one year from the start of the simulation. It does this by subtracting a fixed 365
days — and then correcting:

```r
TmpDateR <- TmpDateR0 - 365 * 24 * 60 * 60
if (format(TmpDateR, format = "%d") != format(TmpDateR0, format = "%d")) {
  TmpDateR <- TmpDateR - 1 * 24 * 60 * 60   # leap year crossed
}
```

If the 365-day step crosses a 29 February, the arrival date slips by one day.
airGR detects this by comparing the day of the month, and removes one more day so
the warm-up starts on the calendar anniversary:

```
2005-01-01  −  365 d  =  2004-01-02     day of month changed
                      →  − 1 d  =  2004-01-01
```

Skip the correction and the warm-up starts one time step late. Every state carried
into the simulation period is then slightly different.

### Why it is easy to miss

The effect is invisible in the obvious test case. A daily simulation starting
1 January 1990 steps back to 1 January 1989 with no leap year in between, so the
correction never fires and a port without it passes.

It surfaced here only on the hourly models, where the demonstration run starts
1 January 2005 and steps back into 2004 — a leap year. GR4H and GR5H were off by
**3 × 10⁻⁶ mm/h**: small, but the single remaining series above tolerance after
everything else had been fixed.

### How GRsuite handles it

[`RunOptions._default_warmup`](../src/grsuite/core.py) reproduces the comparison
on the day of the month, correction included. Passing an explicit `warmup` period
bypasses the whole question.

---

## What both have in common

Neither detail is in a paper, a vignette or a docstring. Both live in the gap
between the model as published and the model as compiled — and both change results
by an amount large enough to fail a strict comparison and small enough to survive
an inattentive one.

That is the argument for validating a port against the reference implementation's
own output, at machine precision, on every commit. A 5 % tolerance would have
passed a port carrying both bugs.

---

## Sources

The two behaviours documented here belong to airGR, not to the published models:

> **Coron, L., Delaigue, O., Thirel, G., Dorchies, D., Perrin, C., Michel, C.**
> airGR: Suite of GR Hydrological Models for Precipitation-Runoff Modelling.
> R package version 1.7.9, INRAE, HYCAR Research Unit, Antony, France.
> doi:[10.15454/EX11NA](https://doi.org/10.15454/EX11NA) —
> files `src/frun_GR4J.f90`, `src/frun_GR6J.f90`, `src/frun_CEMANEIGE.f90` and
> `R/CreateRunOptions.R`.

> **Coron, L., Thirel, G., Delaigue, O., Perrin, C., Andréassian, V.** (2017).
> The suite of lumped GR hydrological models in an R package.
> *Environmental Modelling & Software*, 94, 166–171.
> doi:[10.1016/j.envsoft.2017.05.002](https://doi.org/10.1016/j.envsoft.2017.05.002)

The equations these constants sit inside are published elsewhere — GR4J in Perrin
et al. (2003), GR6J in Pushpalatha et al. (2011), CemaNeige in Valéry (2010).
See **[REFERENCES.md](REFERENCES.md)**.
