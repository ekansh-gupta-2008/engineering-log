# Quarter-Car Ride Model (MATLAB)

2-DOF quarter-car model of a single suspension corner, used to study ride frequency, damping, and the comfort vs. road-holding trade-off. 

## What I built

- A linear 2-DOF model (sprung mass + unsprung mass, spring, damper, tire stiffness) in plain MATLAB. 
- Modal analysis from the eigenvalues of the state matrix (natural frequencies and damping ratios), checked against hand estimates.
- Time simulations with `ode45`: a 5 cm step bump and a sine road at the body natural frequency.
- A damping sweep showing RMS body acceleration (comfort) against RMS tire deflection (road holding).

## Model

States: `x = [zs; vs; zu; vu]` (sprung and unsprung displacement and velocity). Road input `zr`.

```
ms*zs'' = -ks*(zs - zu) - cs*(zs' - zu')
mu*zu'' =  ks*(zs - zu) + cs*(zs' - zu') - kt*(zu - zr)
```

Written as `x' = A*x + B*zr` in `state_matrices()`. Outputs: body acceleration, suspension travel `zs - zu`, tire deflection `zu - zr`.

Wheel-level rates come from the spring and damper via the motion ratio (MR):

```
ks_wheel = k_spring * MR^2
cs_wheel = c_damper * MR^2
```

## Parameters

| Symbol | Value | Meaning | Source |
|---|---|---|---|
| `ms` | 70 kg | Sprung mass, quarter body | Assumed |
| `mu` | 15 kg | Unsprung mass per corner | Assumed |
| `ks` | 18,000 N/m | Suspension stiffness at the wheel | Assumed |
| `cs` | 1,200 N·s/m | Damper coefficient at the wheel | Assumed |
| `kt` | 75,000 N/m | Tire vertical stiffness | Assumed |

Next version replaces `mu` and `ks`/`cs` (via MR) with values taken from the CAD corner.

## Results

**Modal check** (from eigenvalues):

| Mode | Frequency | Damping ratio | Hand estimate |
|---|---|---|---|
| Body | 2.52 Hz | 0.40 | 2.29 Hz |
| Wheel hop | 11.42 Hz | 0.59 | 12.53 Hz |

The hand estimates are about 10% off because they treat the two masses as independent. With `mu/ms` ≈ 0.21 the coupling is not negligible, so the eigenvalues are the reliable numbers.

**5 cm step bump** (`cs` = 1200 N·s/m):
- Body overshoots the 5 cm step by about 46% before settling.
- Peak body acceleration about 30 m/s² (3.1 g).
- Peak suspension travel about 37 mm.
- Final `zs` and `zu` settle at 0.05 m and suspension travel returns to 0, as expected.

**Sine road, 2 cm amplitude at 2.52 Hz:** body amplitude reaches about 1.86 times the road amplitude (resonance amplification).

**Damping sweep** (step input, 3 s window):

| `cs` [N·s/m] | RMS body accel [m/s²] | RMS tire defl. [mm] |
|---|---|---|
| 200 | 3.59 | 5.48 |
| 500 | 2.81 | 3.97 |
| 800 | 2.83 | 3.74 |
| 1200 | 3.10 | 3.86 |
| 2000 | 3.73 | 4.43 |
| 3000 | 4.46 | 5.18 |

Both metrics have a minimum at moderate damping (roughly 500-1200 N·s/m). Too little damping rings after the bump; too much transmits the bump to the body and loads the tire.

Run `quarter_car.m` to reproduce the step response, sine road and damping sweep figures.

## What I learned

- Natural frequency and damping ratio fall out of the eigenvalues of `A`. No need to guess from the time response.
- Simple hand formulas break down as `mu/ms` grows, which matters for hub-motor EVs.

## Limitations

- Linear model: constant spring and damper rates, no bump stops, no asymmetric bump/rebound damping.
- Single step and single sine input. No random road profile yet.
- All parameters are currently assumed, not measured.
