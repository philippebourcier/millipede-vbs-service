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

Code fixes the *level* of $\Delta I$, and phase fixes only how it changes. With float
ambiguities the single-difference ionosphere state reports 30-370 mm, but that figure is
misleading: most of it is a level common to every satellite, which enters every synthesised
phase alike and cancels in the rover's double differences. What reaches a rover is the
*double-differenced* ionosphere. Its sigma is 37-76 mm with float ambiguities, the size of the
1.4-12.5 cm the correction is meant to remove, and **5-9 mm once the cascade fixes them**.
Hence §3.

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

That constraint has a cost, and it is still open. When a later combination's integer follows from
earlier fixes, the constraint fixes it too. For example, the L1-L5 wide lane is the L1-L2 wide
lane plus the L2-L5 extra-wide lane. Its float then sits on an integer before it is tested: such
steps read exactly 0.000 cycles at the moment of rounding. Their ratio test has nothing left to
check, so a wrong earlier fix goes through unchallenged and propagates: one wrong wide-lane cycle
forces a 4-5 cycle narrow-lane error (§9.7). The cascade is therefore only as reliable as its
first fixes, which is why the L1-L2 wide lane, fixed from the UM980's degraded L2 code, matters
so much (§9.8).

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

The floor was raised from 0.5 cm to 2 cm because dividing the whole sigma by $\sqrt n$ claimed a
precision the float does not have: residual troposphere is correlated over minutes, and the
optimism put 1.66 % of integers a cycle out. With the 0.5 cm floor the step did fire. It cut
network misclosure from 0.78 % to 0.01 % at 30 km, yet took the rover from 98.4 % fixed to 13.2 %
(§9.10). Every gain measured afterwards with `--if-narrow-lane` came from what the flag also
does: it **skips** the raw-carrier narrow lane, which cannot work at these distances (§9.10).
Opening the gate instead (to 0.25 cycles) accepts *every* attempted narrow lane, wrong integers
included: 2 244 506 of 2 244 506 below 30 km. Taking the ionosphere from those integers then
loses to the served recipe (3 wins, 21 losses). A phase-levelled ionosphere first needs a
narrow-lane acceptance test that rejects wrong integers, and none exists yet.

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

The interpolated residual itself (the estimated $\Delta\text{ZTD}$ and gradients) is **not
applied** in either served recipe: only the modelled height term is. Applying it took one
Mosaic cluster from 0.56 % to 0.00 % false fixes, but it made mixed and UM980 clusters worse
(on UM980: more harmful distance, 30.9 km against 36.6 km over the ladder). The estimated
residual absorbs whatever else the filter cannot explain.

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

Below about 13 km, direct RTK is already as good as a VBS can make it. Measured on six clean
clusters at 6-13 km over five days, the VBS won 2 cluster-days, lost 5 and traded 23. At 12.77 km
it gave up 88 points of on-spec time (468 against 381, five days summed). Beyond that, whether a VBS helps depends on the family and the distance, and the
failure mode outside the envelope is a confident wrong fix, not a worse float. The rules
(SPEC.md §4) encode the measured envelopes:

- UM980 / UM982: 13-40 km;
- Mosaic-X5: 28-60 km;
- never for all-F9P or mixed clusters.

§7.4 explains why the two families differ.

Nothing measured at run time predicts a wrong fix, so the decision is a prior, taken before a
stream is opened.

### 7.4 Why the two families have different envelopes

| Family | Serve between | Ionosphere correction applied on | GPS L2 used |
|---|---|---|---|
| Unicore UM980 / UM982 | 13-40 km | Galileo, BeiDou | first slot (already L2C) |
| Septentrio Mosaic-X5 | 28-60 km | every constellation | L2C (`--prefer-l2c`) |

Below, "on-spec" is the share of epochs fixed with a de-biased horizontal error of at most
2.5 cm, averaged over five days (§8).

The two receivers are not of different quality. Homogeneous clusters of either family estimate
the per-satellite ionosphere about equally cleanly. Measured as the spread of each satellite's
estimate about its baseline's common level (the level itself is a code bias that cancels at the
rover):

- homogeneous baselines: 2.4-9.4 cm;
- a UM980 differenced against a Mosaic: 21.8-26.2 cm, at only 10-15 km.

That second line is why mixed clusters are refused.

The two bounds have different causes.

**The floor: what the correction adds against what it removes.** The true double-differenced
ionosphere grows with distance: 1.4 cm at 25-40 km, 3.8 cm at 90-110 km (§1.3). The correction's
own error does not shrink with distance. The network fixes wide lanes only (§3.3), so each
satellite's ionosphere state takes its level from code. On some satellites that leaves a
constant error of 15-56 cm, held for hours. A rover's float ambiguity can absorb such a constant;
an integer cannot. Applying the correction pays only where

```math
\underbrace{\nabla\Delta I_\text{true}(d)}_{\text{grows with } d} \;>\; \underbrace{\epsilon_\text{corr}}_{\text{roughly constant}}
```

- **Mosaic.** Its recipe corrects every constellation, including the few satellites with large
  code-levelled errors. Close in, the correction adds more than it removes. Homogeneous Mosaic
  clusters, five days each:

  | Distance to master | On-spec, nearest base | On-spec, VBS | Days won (VBS - nearest) |
  |---|---|---|---|
  | 22.3 km | **70 %** | 53 % | 2 - 2, 1 trade |
  | 26.0 km | **69 %** | 54 % | 1 - 3 |
  | 28.7 km | 53 % | **75 %** | 3 - 2 |
  | 36.9 km | 39 % | **55 %** | 4 - 0 |

  Better base coordinates (§5.3) improve the 22.3 km cluster, but one of its days still had
  wrong fixes on 52 % of fixed epochs, against 1 % for the nearest base. Removing the correction
  below 28 km does not rescue it either. Leave-one-out over four days, the VBS without a
  correction still loses to the nearest base: 80.2 % fixed against 87.4 %, blunders 2.8 % against
  0.8 %, 18 cases won against 42 lost. Close in, a Mosaic nearest base is simply good. Hence
  28 km.
- **UM980.** Its correction is applied on Galileo and BeiDou only. Those estimates are clean
  (correlation with the true ionosphere +0.71 to +0.80), so the correction helps from the start.
  Leave-one-out over four days, it beats the same VBS without a correction in 143 cases against
  81 at 13-40 km, and 46 against 34 at 13-20 km. The floor is then the point below which direct
  RTK is already as good as it gets, 13 km, the same for any family (§7.3).

**The ceiling: how far the correction stays right.** Past 30-40 km the uncorrected ionosphere
alone produces wrong fixes, so a cluster serves only as far as it can correct *every*
constellation the rover uses.

- **Mosaic.** It tracks genuine semi-codeless GPS L2 P(Y), consistent from receiver to receiver,
  and GPS L2C. Its GPS ionosphere is usable, so all constellations are corrected. At 56.9 km a
  Mosaic VBS is on spec 34 % of the time against 9 % for the nearest base, with half the harmful
  distance. Beyond 60 km neither works.
  - The engine reads the first L2 column, which on a Mosaic is often the semi-codeless one.
    Reading its L2C instead (`--prefer-l2c`) is the largest single gain measured on Mosaic
    clusters: harmful distance down 57 %, with 7 % more time on spec, winning 11 of 20
    cluster-days and losing none.
- **UM980.** Its GPS L2 P(Y) code carries a per-unit offset of metres (a 3 m spread across
  units) that drifts by up to metres over days. The GPS ionosphere estimated through it
  correlates poorly with the truth (+0.10 to +0.77), worse than leaving it out (5 cluster-days
  won, 14 lost with GPS added). So GPS satellites, and GLONASS satellites too, get **no**
  ionosphere correction: a third to a half of those in view. Their residual grows with distance
  as in §1.3, and by about 40 km it decides the result:
  - at 38-40 km the UM980 VBS still beats the nearest base (5 cluster-days won, 2 lost, 1 trade);
  - at 51-57 km it only trades, with twice the harmful distance.

  Calibrating the P(Y) offset from each base's own L2C was tried. It does not help: with GPS
  brought to its L2C level, the network ionosphere gets worse, not better.

So the UM980 envelope starts where any VBS can start but stops early, for want of GPS. The
Mosaic envelope reaches further but starts late, because close in neither its corrected nor its
uncorrected VBS beats a Mosaic nearest base.

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
- **Scoring.** A wrong fix and a missing fix do not cost the same: a missing fix is announced
  (an autosteer disengages), a wrong one is silent, and lasts (median 99 s, worst 592 s, that is
  1.6 km of a 10 km/h pass). So every epoch is classed as follows, with the error taken after
  removing each run's median offset:
  - **on-spec**: fixed, horizontal error at most 2.5 cm;
  - **tolerable**: at most 10 cm;
  - **harmful**: more than 10 cm;
  - **no service**: no fix.

  A VBS **beats** the nearest base at a place when, day by day, it gives more on-spec time *and*
  no more harmful distance (distance travelled on harmful epochs). When the two disagree it is a
  **trade**. Days are counted one by one, never pooled first, because pooling five days and then
  testing gives a different answer. A configuration also has to stay under 2 % harmful epochs,
  with no single harmful episode beyond about 120 s. Every envelope was set this way.

---

## 9. Decisions the measurements forced

The choices below each needed real calculation, and each changed the results measurably. They
are grouped by where they act, and every number comes from the engine's test record.

### 9.1 The receiver clock enters the geometry: partial $1 - \dot\rho/c$

An observation tagged at the nominal epoch $t$ was really taken at $t - \delta t_r$, so its range
is $\rho(t - \delta t_r) \approx \rho - \dot\rho\,\delta t_r$. Two free-running receivers drift
apart at about 7.4 µs/s, and one reset itself by exactly 0.9926 ms. A 1 ms step moves the
single-difference geometry by $\dot\rho \times 1$ ms, up to ±0.8 m. It is the same on every band
and different on every satellite, so no slip detector can see it.

The filter already carries $c\,\delta t$ as its clock state. The geometry effect is that same
state with an extra term:

```math
\frac{\partial z}{\partial (c\,\delta t)} = 1 - \frac{\dot\rho}{c}
```

with $\dot\rho$ including the Earth-rotation term. Drift and jumps are then absorbed by the
ordinary update, with nothing to detect and no threshold to tune. On a 29 km baseline over a day,
rejected phase observations fell from 34 653 to 1 359 (39 % → 0.76 %), and observations used rose
from 144 961 to 178 255. Code: `vbs_netekf.c`, `hclk = 1.0 - rdot/CLIGHT`.

### 9.2 The 25° output mask

The corrections for low satellites carry the largest error, and that error grows with distance.
A mask on the *synthesised* stream (not on the network, whose resolution mask stays at 20°)
trades satellites for correction quality. Triple-band clusters deliver about 20 satellites per
VBS epoch and can spare the low ones:

| Case | 20° mask | 25° mask |
|---|---|---|
| 39 km cluster | 96.6 % fixed, 4.68 % false | 98.3 %, 0.30 % |
| 15 km tile corner (10.6 km from the VBS) | 94.9 %, 7.03 % | 98.5 %, 0.52 % |
| 57 km cluster | 52.8 % fixed | 76.0 % |

The 13-fold cut in false fixes at the tile corner is what made a 15 km tile acceptable. Beyond
60 km the mask no longer helps. Code: `out_elmask_deg = 25` in `vbs_session.c`.

### 9.3 Band 6, and a cascade ordered by real wavelength

The combination table is keyed on RINEX band digits, so one row is a different physical
combination on each constellation. Ordered by row, BeiDou ran 0.78 m (B1I-B2a, sitting in an
extra-wide-lane row), then 9.77 m, then 20.93 m (B1C-B1I): shortest first, the wrong way round.
The cascade is now ordered by each combination's actual wavelength,
$\lambda = c/|f_a - f_b|$.

Galileo E6 and BeiDou B3 add extra-wide lanes with real frequency separation: E6-E5b 4.19 m,
B3I-B2b 4.88 m, B2a-B3I 3.26 m, E5a-E6 2.93 m. They also add wide lanes of 0.98-1.02 m. B1C-B1I
is a superb 20.93 m wide lane but only 14 MHz apart: its ionosphere separation,
$\mu_a - \mu_b \approx 0.018$, makes any ionosphere-free use of it ill-conditioned.

On a UM980 cluster, adding band 6 took the VBS from 66.0 % fixed, H95/V95 6.0/8.1 cm and 2.30 %
false to 96.9 %, 4.6/5.4 cm and 0.00 %. Code: `vbs_ambres_cascade_order()`, asserted in
`tests/test_ambres.c`.

### 9.4 Why dual-band RTK is confidently wrong at 20-35 km

A dual-band receiver's longest combination is about 0.86 m (L1-L2 0.862, E1-E5b 0.814, B1I-B2I
0.847; GLONASS adds no extra-wide lane). As the differential ionosphere bias approaches half a
wide-lane cycle, the float is still tight, but centred on the wrong integer. The ratio test
passes, and the fix is confidently wrong. Further out the ratio test fails, and the fixes that
survive are mostly right. The danger band is therefore about 20-35 km: a measured share of all
epochs spent on a wrong fix of 58.9 % at 22.7 km and 28.8 % at 29.1 km, but 0.6 % at 36.9 km.

Extra-wide lanes of 4-21 m keep the bias far from half a cycle, so a triple-band receiver
degrades smoothly instead. At about 30 km of direct RTK:
- dual-band F9P: 86.1 % fixed with 33.46 % false;
- triple-band UM980: 86.3 % with 3.41 % false.

That is the same availability with ten times the wrong fixes. A VBS cut the F9P's wrong epochs to
4.3-5.5 %, against 28.8-37.1 % for direct RTK, which is why F9P clusters were studied rather than
banned outright (§7.3 refuses them for being erratic).

### 9.5 Which L2 column to read

On a Mosaic, the first GPS L2 column is genuine semi-codeless P(Y), and reading its L2C instead
cut harmful distance by 57 % (§7.4). On a UM980 the same switch does nothing: its 2W column
already carries L2C. Measured directly, the UM980's 2W-2X *phase* differs by a median of -9 to
+4 mm, with 8-10 mm scatter on a 244 mm wavelength, while its *code* differs by -0.18 to +0.23 m.
The filter works on phase, so the switch changes which column is read, not which signal.
Code: `ENVELOPE[].prefer_l2c` in `vbs_policy.c`, `vbs_netekf_band_code()`.

### 9.6 Testing whether a correction is right-sized, and why a variance gate cannot

**Paired decomposition.** For a correction $c$ applied to a reference residual, on the same
epochs:

```math
\mathrm{RMS}^2_\text{test} = \mathrm{RMS}^2_\text{ref} - 2\,\mathrm{cov}(\text{ref}, c) + \mathrm{var}(c),
\qquad \alpha = \frac{\mathrm{cov}(\text{ref}, c)}{\mathrm{var}(c)}
```

$\alpha \approx 1$ means the correction is the right size, $\alpha \approx 0$ that it is noise,
and $\alpha < 0$ that it is anti-correlated, which noise cannot produce.
- The dispersive correction scored 0.87 (Galileo) and 0.93 (GPS) at 29 km.
- The estimated non-dispersive residual scored -4.40 and -6.32, with correlation about -0.75.
- The same analysis caught the troposphere model term being counted twice (interpolated, then
  added again).

Withholding the residual repaired a 150 km case: GPS ionosphere-free 6.16 → 4.96 cm, against a
4.92 cm reference. That is why the estimated troposphere residual is not served (§4.2).

**Wrong, or mis-scaled?** A correction off by a constant factor improves steadily as it is
scaled towards the right value. A sweep from 0.85 to 1.50 on a mixed cluster left the vertical
bias flat at -7.7 to -8.0 cm, with a non-monotone fix rate: the correction was wrong, not
mis-scaled. On a UM980 cluster the harm was proportional to how much was applied:

| Scale | Fixed | Blunders |
|---|---|---|
| 1.2 | 55.0 % | 27.67 % |
| 1.0 | 54.7 % | 14.97 % |
| 0.8 | 79.0 % | 0.00 % |
| off | 98.0 % | 0.00 % |

This is where the per-family recipes began.

**A formal variance cannot buy integrity.** Gating corrections on the filter's own variance
(0.05 → 0.025 → 0.012 m) kept the good-to-bad ratio unchanged: about 6 400 correct against 745
wrong became 3 800 against 489. A formal variance measures the filter's confidence, not
systematic error, and a correction poisoned by a wrong integer looks confident.

### 9.7 Reading wrong fixes from their integer signatures

A wide lane wrong by one cycle forces a wrong $N_1 - N_2$, which puts the narrow lane 4-5 L1
cycles off. The data separated perfectly:
- 6 590 narrow-lane errors of that size occurred with the wide lane also misclosed, and none
  without;
- 6 590 is about twice the 3 309 wrong wide-lane fixes;
- 7 486 one-cycle narrow-lane errors occurred with the wide lane closed, and none with it open.

With that propagation removed, the populations have the shape a real wrong-fix population must
have: wide 2 994 : 315 : 0, narrow 7 486 : 26 : 0 (one, two, three or more cycles).

Range and ionosphere errors also leave different per-band signatures. A double-difference range
error $\delta\rho$ misplaces each band by $\delta\rho/\lambda_i$, while an ionosphere error
misplaces it by $\mu_i\,\delta I/\lambda_i$:

| Per metre of error | L1 | L2 | L5 |
|---|---|---|---|
| Range error ($1/\lambda_i$) | 5.25 | 4.10 | 3.93 |
| Ionosphere error ($\mu_i/\lambda_i$) | 5.3 | 6.7 | 7.0 |

Misclosures of 16, 13 and 13 cycles on L1, L2 and L5 matched a 3.1 m range error. The cause was
band-common phase jumps of +55.6 m and then +155.7 m, 30 s apart. Geometry-free is blind to them,
and Melbourne-Wübbena was silent during its warm-up. A fixed 5-cycle Melbourne-Wübbena test
during warm-up, plus resetting a satellite whose phase is rejected 3 epochs running, cut
extra-wide-lane failures from 235 to 73, and raised closure coverage from 78.6 % to 89.2 %.

### 9.8 Measuring ambiguities at the moment of decision

Reading the fractional part $|\hat a - \check a|$ of every float gives +0.000 ± 0.000 over
5.86 million records. That is because an accepted fix is re-applied as a constraint at every
later epoch. Keeping only the first fix of each run leaves 1 762 real decisions. At rounding:
- UM980 GPS L1-L2 wide lane: mean 0.055-0.111 cycles, with 10-27 % of decisions within 0.05 of
  the 0.25 tolerance;
- Mosaic: 0.018-0.034;
- every UM980 combination that avoids band 2: 0.000-0.004.

The median is +0.000 everywhere, so the defect is dispersion, not a bias. It is confined to the
one row built from the L1/L2 Melbourne-Wübbena combination, the only one the UM980's L2 *code*
reaches. Dropping that row for GPS raised delivered usable epochs by 31 % and 49 % on two UM980
geometries, and cost a Mosaic control 10 %. Because the cascade's later steps inherit this row
through the constraint (§3.1), it is the one place where an L2 defect enters every integer.

### 9.9 Tides, and the tide-system convention

Interpolation removes most of a tide, because it is nearly common across a cluster:

| Effect | Raw, at the bases | Residual reaching the rover |
|---|---|---|
| Solid-earth tide, 160 km cluster | 11.7 mm | 0.25 mm |
| Pole tide | 0.34 mm | 0.01 mm |
| Ocean loading, macrotidal coast, 150 km | ~40 mm amplitude | 0.9 mm median, 3.2 mm max |

The ocean-loading residual is predicted by how complex the coastline is, not by the tidal range.

The convention matters more than the effect. The 1005 message announces a conventional
tide-free coordinate, while the synthesis works on the tide-displaced position. Mixing the two
would bias every rover by the permanent tide, 2.4 cm at 43°N to 4.5 cm at 50°N. That bias is
nearly constant in space, so no leave-one-out residual would ever show it. Code:
`vbs_basemeta.c`, `tidedisp(..., TIDE_SOLID_ONLY, ...)`. The library's tide argument is a bit
mask, so its "2" selects ocean loading alone.

### 9.10 The narrow-lane redesign, and why it stopped at the wide lane

The original cascade ended on raw carriers, which needs the slant ionosphere to a fraction of
19-24 cm. It could not hold. On a 46-75 km cluster:
- extra-wide lane: 0 of 59 350 misclosed by a cycle or more;
- wide lane: none beyond 2.5 cycles;
- raw-carrier narrow lane: 22 395 off by one cycle, 4 226 by two, 6 755 beyond 2.5.

The ionosphere-free narrow lane (§3.3) removes the ionosphere from the model, but two things
broke it at the rover:

1. **Reading the ionosphere back from one epoch is too noisy.** The geometry-free read has error
   ```math
   \sigma_{GF} = \frac{2.83\,\sigma_\Phi}{\mu_2 - 1} \approx 2.6\ \text{cm},
   ```
   the size of the signal itself, and correlation with the true ionosphere fell from +0.80 to
   +0.36. Once the integers are fixed and fed back, the filter's own smoothed state *is* the
   geometry-free ionosphere, so no read-back is needed.
2. **The ionosphere-free float depends on the troposphere**,
   ```math
   \hat N = \frac{\nabla\Delta\Phi_{IF} - \Delta m_w\,\mathrm{ZTD} - \text{(known wide-lane term)}}{\lambda_{NL}},
   ```
   so a ZTD error biases every float in an elevation-dependent way that loop closure cannot see
   (it is common to the triangle). Network misclosure fell from 0.78 % to 0.01 % at 30 km, while
   the rover went from 98.4 % fixed to 13.2 %, with 86 % of fixes false and a -12.6 cm vertical
   bias.

That pointed at the troposphere and led to the gradients (§2.2), which fixed the rover. In the
same step the narrow-lane floor rose to 2 cm, and from then on the step never fired (§3.3). The
working VBS therefore resolves wide lanes only.

### 9.11 Precise orbits and Galileo HAS do not help

Satellite clocks and biases cancel between receivers. The orbit error enters only the synthesis,
as

```math
(\mathbf r_\text{SP3} - \mathbf r_\text{brdc}) \cdot (\mathbf e_M - \mathbf e_V) \sim |\delta\mathbf r|\,b/\rho .
```

Measured against final orbits:
- Galileo and BeiDou: under 1 mm at every spacing;
- GPS: 0.76 mm at 25-40 km, 3.83 mm at 140-160 km;
- GLONASS: 2.1-10.5 mm, worst 95 mm, which HAS does not cover.

Removing 4 mm from an 80 mm root-sum-square changes it by 0.1 mm, and ionosphere-free residuals
were identical to 0.1 mm on all five days. Over 280 paired leave-one-out cases the rover gained
+0.12 points of fix rate (95 % CI -0.06 to +0.31), with a sign test of 141/130, p = 0.54. So
`vbsd` takes no orbit-correction stream.

### 9.12 GLONASS: in the network, never in the rover's ambiguity resolution

GLONASS inter-channel biases cancel only between matching hardware (§6.2). Admitting GLONASS to
the network improves every constellation's corrections: blunders fell from 8.86 % to 5.23 %.
But a synthesised stream carries the *master's* hardware signature, not the rover's, and matching
chip models across the bases does not make the rover the same unit. With GLONASS ambiguity
resolution on at the rover, the fixes that survive are 0.00 % wrong against a real base and
63.57 % wrong against a VBS. Rovers should keep GLONASS ambiguity resolution off behind a VBS.

### 9.13 Coordinates must be consistent, not just accurate

§5.3 shows that 7-9 mm of coordinate error changes integers. A partially applied coordinate set
is worse than none: a cluster with some bases moved and others not is internally inconsistent by
however far the two sources disagree. One envelope change was retracted for that reason. Two
reference documents can be compared safely: their disagreement was almost all a common -5.25 mm
vertical shift, which cancels, leaving a median of 1.8 mm horizontal and 3.4 mm vertical (p95
5.0 / 11.1 mm).

Across borders, frames differ in structured ways. One national frame sat at a constant 0.240 m
horizontal translation, flat to 2 mm across latitude. Another's vertical offset grew northward
from 13.8 to 19.9 cm, the signature of post-glacial uplift. Hence the rule that a cluster never
mixes frames.

### 9.14 The statistics of the method

Several measurement traps decided results before any engine change did:

- **Bias versus blunder.** A precise rover sitting 9.2 cm high counted as 32 % "false" under a
  3D > 10 cm rule, on the pass-through arm too. So each run's per-component *median* offset is
  removed before counting blunders. It is the median, because a mean is dragged by the blunders
  it is meant to measure, and per run, never per window, because per-window centring would
  subtract each window's wrong fix. The median breaks down above 50 % blunders, so a split above
  a third is flagged unreliable. On one case, 30.63 % "false" became 13.11 % blunders plus a
  6.4/5.0 cm bias.
- **Choosing test days by the ionosphere, not by Kp.** A Kp 7 "storm" day and a "quiet" day had
  the same double-differenced ionosphere (std 0.038 against 0.039 m). Days are now ranked by the
  robust spread (MAD) of the *change* in double-differenced L1-L2 over a 5-minute lag:
  - the change, not the level, because a cycle slip enters the level as a metre step;
  - MAD, because injecting metre-scale slips on 5 % of epochs left it unchanged;
  - a 5-minute lag, because adjacent 30 s epochs only measure carrier noise, amplified by
    $1/(\mu_2 - 1) = 1.55$.

  The old labels were inverted (storm 0.479, quiet 0.786 cm/min).
- **Denominators.** Raising a code sigma gates more epochs out and flatters the arm under test.
  Counting delivered usable epochs instead showed the default winning by two to three times.
- **Per-day, never pooled.** Pooling five days before a dominance test turned a per-day 3-0-2
  into a "trade". Each day is a test (§8).
- **Confounded labels.** Test clusters were named by rover-to-master distance, but the network
  resolves over master-to-auxiliary baselines, whose mean ran from 61.5 to 164.3 km. A factorial
  grid (4 rover distances × 4 network sizes) separated the two. Read from the engine's own log,
  results track master distance monotonically: 18.7 km 84.7 % on spec, 23.7 km 73.9 %,
  36.5 km 59.6 %.
- **No runtime predictor.** Before a stream opens, distance correlates r = 0.34 with harmful
  share, worst L2 coverage -0.15, coordinate offset -0.06. A crippled station doubles the odds
  of a breach (41 % against 20 %), but refusing on it would wrongly refuse 59 % of cases. The
  rules stay a prior (§7.3).

### 9.15 Auditing the inputs: GPS L2 coverage

A station can declare L2 and fill it on only a few satellites. The useful statistic is GPS
satellite-*records* with L2 over those with L1:
- per epoch, a station sending L2 on 2 of 9 satellites reads 100 %;
- per satellite, the ratio saturates (0.75-0.88 for healthy and broken stations alike);
- the time window must be fixed, because one station read 11 % at 1 Hz and 42 % at 30 s until
  the windows matched.

Of 267 UM980s surveyed, 26 sat below 60 %, 25 of them in an 18.4-22.8 % band. They send L2 only
for the satellites that have no L2C. The cause is firmware: two builds (11833 and 13504) declare
an L2C column and never fill it, while build 17548 is healthy on 77 stations. The defective
stations were inside 161 of 418 test cases, and every envelope that rested on them was
re-measured. `vbs_pipeline` now warns below 60 %: broken stations sit below 23 %, healthy ones
above 76 %.

---

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
