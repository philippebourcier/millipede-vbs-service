# The maths and physics behind the VBS

This document explains what `vbsd` computes and why, from the observation equations to the
bytes a rover receives. It follows the engine's actual models and default constants; where a
measured result decided a choice, the result is given. The numbers come from leave-one-out tests
on a real network of mixed receivers (u-blox F9P, Unicore UM980, Septentrio Mosaic-X5), at
1 Hz and 30 s, on quiet and storm days.

Notation: $s$ a satellite, $i$ a carrier band of frequency $f_i$ and wavelength
$\lambda_i = c/f_i$. $M$ is the master base, $A$ an auxiliary base, $V$ the virtual base. All
quantities are in metres unless stated.

---

## 1. Why RTK degrades with distance

### 1.1 The observation equations

For a receiver $r$ tracking satellite $s$ on band $i$, code $P$ and carrier phase $\Phi$
(in metres) are:

```math
P^s_{r,i} = \rho^s_r + c(\delta t_r - \delta t^s) + T^s_r + \mu_i I^s_r + b_{r,i} - b^s_i + \varepsilon_P
```

```math
\Phi^s_{r,i} = \rho^s_r + c(\delta t_r - \delta t^s) + T^s_r - \mu_i I^s_r + \lambda_i N^s_{r,i}
+ \beta_{r,i} - \beta^s_i + \varepsilon_\Phi
```

- $\rho$ is the geometric range, including the Sagnac correction;
- $\delta t_r$, $\delta t^s$ are the receiver and satellite clocks;
- $T$ is the troposphere (non-dispersive: the same on every band);
- $I$ is the slant ionosphere on a reference frequency $f_1$, scaled to band $i$ by
  $\mu_i = (f_1/f_i)^2$. First-order ionosphere delays the code and advances the phase, so the
  sign flips;
- $N$ is the integer ambiguity;
- $b$ and $\beta$ are code and phase hardware biases;
- code noise $\varepsilon_P$ is decimetres, phase noise $\varepsilon_\Phi$ is millimetres.

### 1.2 Double differences

RTK differences between two receivers (rover $R$, base $B$) and then between two satellites $s$
and a reference $q$:

```math
\nabla\Delta(\cdot) = \left[(\cdot)^s_R - (\cdot)^s_B\right] - \left[(\cdot)^q_R - (\cdot)^q_B\right]
```

Clocks and hardware biases cancel, and the ambiguity $\nabla\Delta N$ is an integer:

```math
\nabla\Delta\Phi_i = \nabla\Delta\rho + \nabla\Delta T - \mu_i \nabla\Delta I
+ \lambda_i \nabla\Delta N_i + \varepsilon
```

What does not cancel is the **spatial decorrelation** of the atmosphere and of the orbit error
over the baseline: $\nabla\Delta I$, $\nabla\Delta T$ and a term of about $|\delta \mathbf{r}^s|\,
d/\rho$ for an orbit error $\delta\mathbf{r}^s$ over a baseline $d$. Broadcast orbits are good to
about a metre, worth 2.5 mm at 50 km: small. The atmosphere is not.

### 1.3 How big the residual is

Measured as double-difference residuals between real bases, uncorrected, after removing each
arc's constant (quiet day, 15° mask):

| Baseline | Geometry-free (ionosphere) | Ionosphere-free (troposphere + orbit) |
|---|---|---|
| 25-40 km | 1.4 cm | 8.2 cm |
| 50-70 km | 2.4 cm | 5.9 cm |
| 90-110 km | 3.8 cm | 8.7 cm |
| 140-160 km | 5.1 cm | 7.6 cm |

On a geomagnetic storm day (Kp 8.7) the ionospheric column grows by 2.1-2.4 times, to 12.5 cm at
150 km. An L1 cycle is 19 cm. A rover deciding integers at a quarter of a cycle has very little
margin past 20-30 km. Worse, it can decide them wrongly while reporting a centimetre sigma.

The VBS idea: estimate these residuals on the network, where the base coordinates are known,
and remove them for the rover.

---

## 2. The network filter (one per master-auxiliary baseline)

`vbs_netekf` runs an extended Kalman filter per baseline $M \to A$. Both ends are surveyed
points, so there is **no position state**. That is what keeps a few-hundred-state filter well
conditioned over 150 km.

### 2.1 Single differences, not double

The filter works on **between-station single differences** ($A - M$, per satellite), not double
differences. In a double-difference state every quantity refers to one reference satellite.
When that satellite sets, the whole state must be re-expressed, and any error in that transform
becomes a step in the corrections. Single differences have no reference satellite, so the
corrections are continuous by construction. The cost is one receiver-clock-difference state per
constellation, and single-difference ambiguities that are not integers (they absorb the two
receivers' phase biases). Integer structure is only needed at resolution time, where the double
differences are formed and the biases cancel.

After removing the computed range difference $\Delta\rho$ and an a priori troposphere
difference $\Delta T_\text{model}$ (Saastamoinen, 70 % humidity), per satellite and band:

```math
z_P = \Delta P - \Delta\rho - \Delta T_\text{model}
= \Delta\text{clk} + m_w\,\Delta\text{ZTD} + \mathbf{m}_g\cdot\Delta\mathbf{G} + \mu_i\,\Delta I + \varepsilon_P
```

```math
z_\Phi = \Delta\Phi - \Delta\rho - \Delta T_\text{model}
= \Delta\text{clk} + m_w\,\Delta\text{ZTD} + \mathbf{m}_g\cdot\Delta\mathbf{G} - \mu_i\,\Delta I
+ \lambda_i\,\Delta N_i + \varepsilon_\Phi
```

The opposite sign of $\Delta I$ on code and phase separates the ionosphere from the ambiguity.
The spread of $\mu_i$ across bands separates it from the non-dispersive terms.

### 2.2 State vector

| State | Count | Process noise (default) |
|---|---|---|
| receiver clock difference $\Delta\text{clk}$ | 1 per constellation | re-seeded every epoch from the median code single difference |
| zenith troposphere difference $\Delta\text{ZTD}$ | 1 | random walk $10^{-8}$ m²/s (≈6 mm/h) |
| horizontal gradients $\Delta G_N, \Delta G_E$ | 2 | random walk $10^{-12}$ m²/s, initial σ 2 mm |
| slant ionosphere difference $\Delta I^s$ | 1 per satellite | random walk $10^{-4}$ m²/s as served (library default $10^{-6}$) |
| single-difference ambiguity $\Delta N^s_i$ | 1 per satellite and band (up to 5 bands) | constant; reset on a slip |

Wet mapping $m_w$ is the Niell function. The gradients use the Chen-Herring mapping

```math
m_g(e) = \frac{1}{\sin e \tan e + 0.0032}, \qquad
\mathbf{m}_g\cdot\Delta\mathbf{G} = m_g(e)\,(\Delta G_N\cos a + \Delta G_E \sin a)
```

capped at 10 so that the lowest satellites cannot dominate.

**Why gradients matter.** A single zenith scalar cannot represent a tilted atmosphere. What it
leaves behind depends on azimuth and elevation, is worth centimetres on low satellites, and
biases every baseline of a cluster alike. Loop closure never sees it, but the corrections carry
it. Left loose, gradients are catastrophic. Tightly constrained, they cut false fixes threefold.

### 2.3 Weighting and fault exclusion

Zenith sigmas are 0.6 m for code and 6 mm for phase, with variance $\sigma^2(1 + 1/\sin^2 e)$.
Observations are
processed one at a time. Each is tested against the predicted innovation variance before it
touches the state, and rejected beyond 4σ (after 10 epochs of convergence). A satellite whose
phase is rejected for 3 epochs in a row has its ambiguities reset. A fix the phase no longer
supports is not a fix.

What the filter deliberately does **not** estimate is a per-satellite non-dispersive term:
orbit error and any asymmetric troposphere go into the ambiguities. At these baseline lengths
that is a few millimetres.

### 2.4 Why the ambiguities must be fixed

With float ambiguities, the slant ionosphere converges only to **30-370 mm**. Code fixes the
*level* of $\Delta I$, and phase fixes only how it changes. That is an order of magnitude above
the 1.4-12.5 cm the correction is meant to remove. Fixing the ambiguities turns $\Delta I$ into
a phase-determined quantity. Hence §3.

---

## 3. Network ambiguity resolution

### 3.1 The cascade

Direct L1/L2 fixing at 50-150 km is not realistic: an L1 ambiguity is 19 cm, so the float must
be good to a few centimetres. The way through is to fix long-wavelength combinations first. For
two bands $a, b$:

```math
\lambda_{WL} = \frac{c}{f_a - f_b}, \qquad
\Phi_{WL} = \frac{f_a\Phi_a - f_b\Phi_b}{f_a - f_b}
```

| Step | Example combinations | Wavelength |
|---|---|---|
| extra-wide lane | E5a-E5b, B2a-B2b | ≈ 9.8 m |
| wide lane | L1-L2, E1-E5b, L1-L5 | ≈ 0.75-0.86 m |
| narrow lane | the individual bands | ≈ 0.19-0.25 m |

The steps are ordered by each combination's actual wavelength on each constellation. Galileo
E6 and BeiDou B3 are used where the bases stream them. Each fixed step is applied back to the
filter as a tight constraint (variance $10^{-6}$ cycles²), so the next step starts from a float
the previous one already sharpened. Triple-frequency is the normal case between Unicore or
Septentrio bases. F9P baselines are dual-frequency and start at the wide lane.

Over the integer candidates, the search is LAMBDA (MLAMBDA). Validation is the ratio test:
second-best over best squared norm, at least 3. When the test fails, the worst-conditioned
ambiguity is dropped and the subset retried (partial resolution, minimum 2 fixed). Satellites
below 20° are left out of resolution.

### 3.2 Time correlation at 1 Hz

The filter's covariance assumes independent epochs. At 1 Hz consecutive residuals are
strongly correlated (multipath, troposphere), so the formal ambiguity covariance is too
optimistic by roughly the correlation time over the interval. The search, and only the search,
sees $Q$ inflated by $\tau/\Delta t = 30\text{ s}/1\text{ s}$. The float keeps converging at
full rate, while the ratio test sees the information the data really holds.

### 3.3 The ionosphere-free narrow lane

With the wide-lane integer $N_{WL} = N_a - N_b$ known, the remaining integer can be resolved on
the ionosphere-free combination, which removes the ionosphere and its baseline-length limit:

```math
\Phi_{IF} = \frac{f_a^2\Phi_a - f_b^2\Phi_b}{f_a^2 - f_b^2}
= \rho + T + \lambda_{NL} N_a + \frac{c f_b}{f_a^2 - f_b^2} N_{WL}, \qquad
\lambda_{NL} = \frac{c}{f_a + f_b}
```

One epoch of $\Phi_{IF}$ is about a third of a narrow-lane cycle, so its float is averaged over
the arc (at least 120 epochs). It is accepted only when:

- its sigma is below 0.08 cycles, including a 2 cm systematic floor (troposphere residual)
  that averaging cannot remove;
- its mean is within 0.12 cycles of an integer.

**As shipped, this step does not fire.** The 2 cm floor over a 10.7 cm narrow lane is
0.19 cycles, above the 0.08 gate. In practice the network fixes **wide lanes only**, and each
auxiliary's ionosphere takes its level from code. That has consequences (§5.3).

### 3.4 Loop closure

In a cluster, $M\to A_1$, $M\to A_2$ and $A_1\to A_2$ must agree: their fixed integers must close
to zero around every triangle. The resolver reports the worst misclosure, and the offline
pipeline (`tools/vbs_pipeline`) checks it. A misclosure flags a wrong fix, or a slip that one
filter saw and another missed. That second failure is why slips are detected once per *station*
(§6), never per baseline. Closure is blind to a per-satellite constant common to all baselines,
such as an unmodelled troposphere gradient (§2.2).

---

## 4. From baselines to the virtual point: interpolation

Each auxiliary $k$ delivers, per satellite, a correction relative to the master,
$\Delta I^s_k$ and $\Delta T^s_k$, with its variance. These are samples of a field, known at
the auxiliaries and zero at the master by definition. The VBS needs the field at $V$.

### 4.1 Least-squares collocation, anchored at the master

The classic method is Wanninger's weighted plane, $\Delta(E,N) = a + bE + cN$. Its residual
grows by about 0.19 cm per km of master-to-VBS displacement, because a plane cannot carry
curvature. The engine uses **least-squares collocation** instead. With samples $\mathbf{s}$
(the auxiliaries, plus the master as a sample of value zero):

```math
\hat{s}_V = \mathbf{c}_0^\top (\mathbf{C} + \mathbf{N})^{-1}\mathbf{s}, \qquad
\sigma^2_V = \sigma_f^2 - \mathbf{c}_0^\top (\mathbf{C} + \mathbf{N})^{-1}\mathbf{c}_0
```

with an isotropic exponential covariance

```math
C(d) = \sigma_f^2\, e^{-d/L}, \qquad \sigma_f = 5\text{ cm}, \quad L = 60\text{ km}
```

$\mathbf{C}$ holds the covariance between samples, $\mathbf{c}_0$ the covariance between each
sample and $V$, and $\mathbf{N}$ the samples' own noise (diagonal, from the filter).
Collocation carries curvature without more parameters than samples. It also weights a distant
or noisy auxiliary by how much the field at $V$ is expected to resemble it. On a 10 km tile
corner it took a rover from 90.5 % fixed with 5.9 % false to 98.0 % with 0.6 %.

**The master anchor** matters on its own. A free fit through three auxiliaries discards the
known zero at the master and lets an offset and tilt through, which reaches the rover as a
horizontal bias growing with the displacement.

A sample with σ above 5 cm is ignored. A satellite needs three primary auxiliaries, or two
plus a **reserve**: one more station of the same family within 250 km, fitted for a satellite
only when the primaries fall short ("2 + 1 backup"; `--no-backup` requires three).

### 4.2 The height term is modelled, not interpolated

Troposphere depends strongly on height. The a priori model difference between master and VBS
heights is computed (Saastamoinen at a common elevation), and only the *residual* is
interpolated. Interpolating the whole term across stations at different altitudes would smear a
height signal across a horizontal plane.

### 4.3 Extrapolation

The leverage, the variance inflation of the prediction at $V$ from geometry alone, is reported
with every correction. Inside a well-spread cluster it is about $1/n$ (0.25-0.35 measured); just
outside the hull it reaches 6-12. Prediction degrades quickly outside the cluster.

---

## 5. Synthesising the virtual base

### 5.1 The translation

The VBS is the master's own observations moved to $V$:

```math
\Delta\rho^s = |\mathbf{X}^s - \mathbf{X}_V| - |\mathbf{X}^s - \mathbf{X}_M|
```

```math
\Phi^s_{V,i} = \Phi^s_{M,i} + \Delta\rho^s + \Delta T^s - \mu_i\,\Delta I^s, \qquad
P^s_{V,i} = P^s_{M,i} + \Delta\rho^s + \Delta T^s + \mu_i\,\Delta I^s
```

with $\Delta I^s$, $\Delta T^s$ the interpolated corrections at $V$, and phase converted to
cycles by each signal's own wavelength.

- **The integer nature is preserved exactly.** Every term added is a known quantity in metres.
  Nothing introduces a fractional bias, so the rover's double-difference ambiguity against the
  VBS is an integer.
- **The light-time is redone for $V$.** The signal reaching $V$ at the same reception epoch was
  sent about 100 µs earlier or later over 30 km, during which the satellite moves about 0.4 m.
  The satellite position is propagated to $V$'s own transmission time and the range
  re-converged, with the Sagnac correction on both ranges.
- A satellite with no valid correction is not output at all.

### 5.2 What the rover then sees

Against the VBS, the rover's double differences are those of a short baseline. The residual is
the interpolation error, not the full atmosphere over 30-60 km. The rover is unchanged: any RTK
receiver accepts the stream as an ordinary base.

### 5.3 What limits it in practice

**Receiver code biases.** A single-difference ionosphere state cannot separate the true
differential ionosphere from the two receivers' differential code biases: with wide lanes only,
its level comes from code. Measured on real clusters, this per-pair constant is 0.5-2.5 m. It
differs per auxiliary and per constellation, and drifts 16-95 cm over a day. Interpolating
three such constants mostly interpolates receiver biases. The datum can be removed per
auxiliary (`--dd-datum`). That improves the dispersive correction, but on UM980 clusters it
still lost to serving no ionosphere at all. The next point is what finally worked there.

**The UM980's GPS L2.** Its semi-codeless P(Y) code carries a per-unit offset of metres from
its L2C (3 m spread across units), drifting over days. The GPS network ionosphere on UM980
clusters therefore correlates poorly with the truth (+0.10 to +0.77), while Galileo's stays at
+0.71 to +0.80. The served recipe applies the dispersive correction on **Galileo and BeiDou
only** for UM980 clusters. Calibrating the P(Y) offset from the base's own L2C was tried and
does not help.

**Families.** Hardware biases cancel in a double difference only between like receivers. A
cluster mixing families was never better than the nearest base. All-F9P clusters (dual-band,
no third frequency) are erratic: they sometimes win, and nothing predicts when.

**Coordinates.** A coordinate error at a base enters every satellite's correction as a
satellite-dependent geometric error. Moving base coordinates by 7-9 mm, to a precise
point-positioning solution, **halved** a rover's distance travelled on wrong fixes, even on
clusters where no base was off by more than 9 mm. Seven millimetres is 4 % of an L1 cycle, and
five stations' worth feeds a decision taken at 0.25 cycles. Hence the coordinate document
(SPEC.md §5.2), and the gate refusing a base whose declared position sits more than 10 mm
horizontally or 20 mm vertically from where its neighbours put it.

**Earth tides.** Solid-earth tide (up to ~30 cm vertically) is applied to every base
coordinate. Across a 50 km cluster it is almost common, but not exactly.

---

## 6. Data hygiene: slips, signals, lock time

### 6.1 Cycle slips, per station

Three detectors, each blind to something the others catch:

- **Geometry-free**, $L_a - L_b$ in metres: cancels geometry, clocks and troposphere, leaving
  ionosphere and the ambiguity difference. The threshold grows with the interval,
  $0.03\text{ m} + 0.004\text{ m/s}\cdot\Delta t$, because the ionosphere moves more between
  30 s epochs than between 1 s epochs. GF is blind to slips like 9:7 cycles on L1:L2, which
  shift it by under a millimetre.
- **Melbourne-Wübbena**, wide-lane phase minus narrow-lane code,
  $MW = \frac{f_a\Phi_a - f_b\Phi_b}{f_a - f_b} - \frac{f_a P_a + f_b P_b}{f_a + f_b} = \lambda_{WL} N_{WL} + \varepsilon$.
  It cancels geometry, clocks, troposphere and ionosphere, leaving the wide-lane ambiguity. It is
  noisy (code), so it is tested against a running mean: max(1.5 cycles, 4σ). A fixed 5-cycle test
  bridges the warm-up after each slip.
- **Lock time / loss-of-lock** from the receiver. It is kept as a trigger but never trusted
  alone: on a storm day it reported zero slips where GF and MW found 1 991.

A gap of more than 3 intervals ends the arc. The arc id per station and satellite is the single
source of truth for every filter using that station.

### 6.2 Signals matched by band, not by label

Two receivers label the same physical signal differently (GPS L2C as 2X, 2L or 2S). Signals are
matched by (constellation, RINEX band digit), never by tracking attribute. GLONASS (FDMA) is
handled per satellite frequency, $f = f_0 + k\,\Delta f$ for channel $k$. Its inter-channel
biases cancel in a double difference only when both ends have matching hardware delays, so
GLONASS enters the network only on homogeneous clusters, with every receiver of the same
family. Otherwise it is excluded.

### 6.3 Telling the rover about discontinuities

MSM has no loss-of-lock flag; the **lock-time indicator** is the only signal. The VBS resets it
when:

1. the master slips;
2. a network ambiguity is re-fixed or lost, which steps the correction by whole cycles;
3. the interpolated correction jumps, for instance when an auxiliary drops out;
4. a satellite returns after exclusion.

Cases 2-4 are why this is not a pass-through of the master's lock time: the VBS carries a
*computed* phase. A missed reset makes the rover keep a fix that is wrong by whole cycles while
its sigma stays confident.

---

## 7. Serving geometry

### 7.1 One synthesis per tile

The VBS point is the centre of a fixed global tile. Every rover on a tile shares one synthesis,
so cost scales with occupied tiles, not with users, and the worst rover-to-VBS distance is what
the tile shape trades against the tile count.

A hexagon is the most compact cell that tiles a plane. With circumradius $R$, it covers
$\tfrac{3\sqrt3}{2}R^2 \approx 2.60R^2$, while a square with the same worst case (half-diagonal
$R$) covers $2R^2$. At $R = 10.6$ km that is 292 km² against 225 km²: the same worst-case distance
with about 23 % fewer tiles. Only three hexagons meet at a vertex (four squares do), and the six
neighbours are equidistant.

On the sphere the hexagonal grid is built as rows of equal latitude, $1.5R$ apart, with centres
$\sqrt3 R$ apart along a row and odd rows shifted by half a column. A position belongs to the
nearest centre, so the cells are the Voronoi cells of that lattice. The column count must shrink
with $\cos\varphi$ towards the poles. If each row had its own count (about 4 columns fewer per
row at 48°), the half-column offset between rows would drift around the globe, and wherever it
reached zero the cells would turn into rectangles. Counts are therefore held constant over bands
of rows while a cell's width stays within 5 % of nominal. Only the seam between two bands is
irregular. Measured:

- 97-98 % of cells have six neighbours;
- 99 % of positions are within 10.0 km of their VBS, and the worst case is 11.4-11.7 km, on seam
  rows only;
- the grid needs 21 % fewer tiles than 15 km squares.

The tile's RTCM station ID differs from every neighbour's: a rover crossing a boundary sees its
base change, which is how RTCM signals a new reference. A 300 m hysteresis, measured past the
bisector with the new cell, prevents flapping on a boundary.

### 7.2 Pass-through while converging

A new cluster needs about three minutes before its first correction passes its gate. In the
meantime the rover gets the nearest base's own stream, verbatim, under that base's own position.
A translation without corrections would carry the master-to-tile atmosphere under a
short-baseline position, and a rover weighting its ionosphere by baseline length would fix wrong
integers on it.

### 7.3 When to serve at all

Direct RTK is already at the acceptance ceiling below about 13 km: 97.7 % fixed with 1.2 % false
at 10.7 km. Beyond that, whether a VBS helps depends on the family and the distance, and the
failure mode outside the envelope is a confident wrong fix, not a worse float. The rules
(SPEC.md §4) encode the measured envelopes:

- UM980 / UM982: 13-40 km;
- Mosaic-X5: 28-60 km;
- never for all-F9P or mixed clusters.

Nothing measured at run time predicts a wrong fix, so the decision is a prior, taken before a
stream is opened.

---

## 8. How the claims were tested

- **Leave-one-out.** Remove a real base, synthesise a VBS at its exact coordinates from its
  neighbours, and double-difference the result against the base's own observations. Geometry
  cancels exactly, so the residual is the correction error (`tools/vbs_loo`). The same
  synthesised stream is then fed to a rover engine, and its fix rate, false fixes and
  "harmful distance" (distance travelled on a wrong fix) are compared with the nearest base
  under the same conditions.
- **Zero baseline.** Synthesising with zero corrections at the master's own position must
  reproduce the master's observations exactly. Through the MSM7 encoder and an independent
  decoder (RTKLIB), they must match within one MSM least significant bit: 0.56 mm on code,
  0.14 mm on phase (`tests/test_zero_baseline.c`).
- **Acceptance**, relative to the alternative at the same place: at least 90 % fixed, H95/V95
  at most 10 cm, at most 1 % false fixes.

## References

- Wanninger, L. (1995). Improved ambiguity resolution by regional differential modelling of the
  ionosphere. *Proc. ION GPS-95*.
- Vollath, U., Buecherl, A., Landau, H., Pagels, C., Wagner, B. (2000). Multi-base RTK
  positioning using virtual reference stations. *Proc. ION GPS 2000*.
- Teunissen, P. J. G. (1995). The least-squares ambiguity decorrelation adjustment: a method
  for fast GPS integer ambiguity estimation. *Journal of Geodesy* 70, 65-82.
- Chang, X.-W., Yang, X., Zhou, T. (2005). MLAMBDA: a modified LAMBDA method for integer
  least-squares estimation. *Journal of Geodesy* 79, 552-565.
- Moritz, H. (1980). *Advanced Physical Geodesy* (least-squares collocation).
- Chen, G., Herring, T. A. (1997). Effects of atmospheric azimuthal asymmetry on the analysis of
  space geodetic data. *JGR* 102(B9).
- Niell, A. E. (1996). Global mapping functions for the atmosphere delay at radio wavelengths.
  *JGR* 101(B2).
- Saastamoinen, J. (1972). Atmospheric correction for the troposphere and stratosphere in radio
  ranging of satellites. *Geophysical Monograph* 15.
- Melbourne, W. G. (1985); Wübbena, G. (1985). The wide-lane / narrow-lane code-phase
  combination, *Proc. 1st Int. Symp. Precise Positioning with GPS*.
- Petit, G., Luzum, B. (eds.) (2010). *IERS Conventions*, ch. 7 (solid-earth tide).
