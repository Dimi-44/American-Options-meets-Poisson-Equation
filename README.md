# Pricing American Options through the Poisson Equation

Code for my BSc project at the Athens University of Economics and Business (2024), *"Applying the Poisson equation to the pricing of American derivatives"*. The full write-up, in Greek, if you want to access it feel free to contact me. 

The idea of the project is that pricing an American option can be recast as a problem with the same structure as the Poisson equation. Once it is in that form, the numerical tools developed for Poisson apply directly: finite differences to discretise it and Successive Over-Relaxation (SOR) to solve the resulting linear systems.

## 1. The problem

An option gives its holder the right, but not the obligation, to buy (call) or sell (put) an asset at a fixed strike price $K$.

* A **European** option can only be exercised at maturity $T$. Under the Black–Scholes assumptions its value $V(S,t)$ solves

$$
\frac{\partial V}{\partial t} + \tfrac12\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + rS\frac{\partial V}{\partial S} - rV = 0,
$$

  and has a closed-form solution.

* An **American** option can be exercised at **any time** up to $T$. The holder exercises as soon as the payoff $\Phi(S)$ is worth more than keeping the option alive, where $\Phi(S) = \max(K-S,0)$ for a put and $\max(S-K,0)$ for a call. The value must therefore always satisfy $V \ge \Phi$.

This early-exercise right turns pricing into a **free-boundary problem**. For every time $t$ there is an unknown critical price $S^\ast(t)$ that separates two regions:

* in the **continuation region**, the option is worth more alive and the Black–Scholes equation holds;
* in the **exercise region**, $V = \Phi$.

The boundary $S^\ast(t)$ is part of the solution, so no closed-form price exists and numerical methods are needed.

## 2. The Poisson connection

### The Poisson equation and the obstacle problem

The Poisson equation $-\nabla^2 u = f$ describes, for example, the shape $u$ of an elastic membrane (a sheet, a drum skin) pushed by a load $f$ and fixed at its edges.

Now place an **obstacle** $\psi$ under the membrane, so that it can rest on it but not pass through it. This is the classical **obstacle problem**:

$$
-\nabla^2 u \ge f, \qquad u \ge \psi, \qquad \bigl(-\nabla^2 u - f\bigr)\,(u - \psi) = 0 .
$$

At every point, either the membrane is free and Poisson holds, or it touches the obstacle and $u = \psi$. Where the contact region ends is not known in advance. It is a free boundary.

### The American option has the same structure

Write the Black–Scholes operator in time-to-maturity $\tau = T - t$:

$$
\mathcal{L}V = \frac{\partial V}{\partial \tau} - \tfrac12\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} - rS\frac{\partial V}{\partial S} + rV .
$$

The American option is then

$$
\mathcal{L}V \ge 0, \qquad V \ge \Phi, \qquad \mathcal{L}V\,(V - \Phi) = 0 ,
$$

which is the obstacle problem again, with the payoff playing the role of the obstacle. Black–Scholes is an evolution (parabolic) equation, so the match is exact once time is discretised. Each implicit time step

$$
\bigl(I + \Delta t\,\mathcal{A}\bigr) V^{n} \ge V^{n-1}, \qquad V^{n} \ge \Phi, \qquad \text{(complementarity)}
$$

is an obstacle problem for a Poisson-type (elliptic) operator. The table shows the correspondence.

| Poisson obstacle problem | American option (one time step) |
|---|---|
| membrane height $u(x)$ | option value $V(S)$ |
| operator $-\nabla^2$ | $I + \Delta t\,\mathcal{A}$, with $\mathcal{A} = -\tfrac12\sigma^2S^2\partial_{SS} - rS\,\partial_S + r$ |
| load $f$ | option values from the previous step, $V^{n-1}$ |
| obstacle $\psi$ | payoff $\Phi(S)$ |
| contact region $u = \psi$ | exercise region $V = \Phi$ |
| edge of the contact region | early-exercise boundary $S^\ast$ |
| fixed edges of the membrane | values at $S = 0$ and $S = S_{\max}$ |

Because the two problems share this structure, the method used for the obstacle problem, **projected SOR**, prices the American option.

## 3. Numerical scheme

**Grid.** Prices $S_i = i\,\Delta S$ for $i = 0,\dots,M$ on $[0, S_{\max}]$, and time steps $\tau_n = n\,\Delta t$ for $n = 0,\dots,N$.

**Finite differences.** Use backward Euler in time and central differences in $S$. Since $S_i/\Delta S = i$, each time step gives a **tridiagonal** system, just like the discrete Poisson equation:

$$
a_i V_{i-1}^{n} + b_i V_i^{n} + c_i V_{i+1}^{n} = V_i^{n-1},
\qquad
\begin{aligned}
a_i &= -\tfrac{\Delta t}{2}\,(\sigma^2 i^2 - r i)\\
b_i &= 1 + \Delta t\,(\sigma^2 i^2 + r)\\
c_i &= -\tfrac{\Delta t}{2}\,(\sigma^2 i^2 + r i)
\end{aligned}
$$

**Projected SOR.** Sweep through the interior nodes. At each node, take a Gauss–Seidel step, over-relax it by $\omega \in (1,2)$, and project the result back above the payoff:

$$
V_i \leftarrow \max\!\Bigl( (1-\omega)\,V_i + \frac{\omega}{b_i}\bigl(V_i^{n-1} - a_i V_{i-1} - c_i V_{i+1}\bigr),\; \Phi(S_i) \Bigr)
$$

Sweeps repeat until the largest change in a sweep falls below a tolerance. The projection, the `max(…, payoff)`, is what enforces early exercise; without it the scheme would price a European option.

**Boundary conditions.**

| | $S = 0$ | $S = S_{\max}$ |
|---|---|---|
| Put | $V = Ke^{-r\tau}$ | $V = 0$ |
| Call | $V = 0$ | $V = S_{\max} - Ke^{-r\tau}$ |

## 4. Code

The same algorithm is implemented in both languages:

* [`R/american_sor.R`](R/american_sor.R), compared with **RQuantLib**
* [`python/american_sor.py`](python/american_sor.py), compared with **QuantLib**

Each script prices an at-the-money put and call ($S_0 = K = 100$, $T = 1$, $r = 5\%$, $\sigma = 20\%$) with two settings:

* **Base SOR:** $M = 100$, $N = 1000$, $\omega = 1.2$.
* **Optimised SOR:** $M = 300$, $N = 3000$. The grid is aligned so the strike sits exactly on a node, and $\omega = 2/(1+\sin(\pi/M))$, the optimal relaxation parameter for the discrete Poisson problem.

### Running

```bash
# Python
pip install -r python/requirements.txt
python python/american_sor.py

# R (installs RQuantLib if missing)
Rscript R/american_sor.R
```

The loops are written out explicitly for readability, so the optimised run is slow: about 19 minutes in Python and 12 in R.

### Output (Python)

```
[PUT OPTION]
Base SOR Price (M=100, N=1000):  6.0984  | Diff vs QL: 0.011386
Opt SOR Price (M=300, N=3000):   6.0870  | Diff vs QL: 0.000001
QuantLib Reference Price:       6.0870

[CALL OPTION]
Base SOR Price (M=100, N=1000):  10.4647  | Diff vs QL: 0.012504
Opt SOR Price (M=300, N=3000):   10.4478  | Diff vs QL: 0.004395
QuantLib Reference Price:       10.4522
```

### Output (R)

```
Type   | Base SOR     | Opt SOR      | RQuantLib    | Err (Base)   | Err (Opt)
------------------------------------------------------------------------------
Put    | 6.0984       | 6.0870       | 6.0976       | 0.000786     | 0.010601
Call   | 10.4647      | 10.4478      | 10.4506      | 0.014078     | 0.002821
```

Both languages give identical SOR prices.

## 5. Notes from a later review

* **Accuracy.** More accurate reference prices are 6.0904 for the put (high-resolution finite differences) and 10.4506 for the call. The call value is exact: with no dividends an American call is never exercised early, so it equals the Black–Scholes call. Against these references the base grid is off by about 0.01 and the refined grid by about 0.003. The QuantLib reference in the Python script is itself a coarse 200×200 grid, which is why the optimised put looks almost exact.
* **RQuantLib.** `AmericanOption()` uses the Barone-Adesi–Whaley approximation by default and ignores `timeSteps`. That approximation overprices the put by about 0.007, which is why the R table makes the base put look more accurate than the optimised one; the opposite is true. Pass `engine = "CrankNicolson"` for a finite-difference price.
* **Where the Poisson analogy stops.** The Black–Scholes matrix is far more diagonally dominant than the Poisson matrix, so the Poisson-optimal $\omega \approx 1.98$ does not carry over. The optimal value here is about 1.08, which gives the same prices with roughly 85 times fewer SOR sweeps.
* **Put boundary at $S = 0$.** For an American put the exact value there is $K$, since the holder would exercise immediately. The code uses the European value $Ke^{-r\tau}$, which has no visible effect at $S_0 = 100$.

## References

* Seydel, R. (2006). *Tools for Computational Finance* (3rd ed., Ch. 4). Springer.
* Hull, J. C. (2018). *Options, Futures, and Other Derivatives* (10th ed.). Pearson.
* Wilmott, P., Howison, S., & Dewynne, J. (1995). *The Mathematics of Financial Derivatives*. Cambridge University Press.
* Cryer, C. W. (1971). The solution of a quadratic programming problem using systematic overrelaxation. *SIAM Journal on Control*, 9(3), 385–392.
* Pinsky, M. A. (2003). *Partial Differential Equations and Boundary-Value Problems with Applications* (3rd ed.). AMS.

## License

MIT, see [LICENSE](LICENSE).
