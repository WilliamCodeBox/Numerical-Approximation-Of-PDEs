# Numerical Approximation Of PDEs

A project for numerical approximation of Partial Differential Equations with C++

A good starter project for entering the world of scientific computing and numerical computing.

PDEs are mathematical models of real world physical problems. Following three kinds of PDE are the most classical and basic models.

- Convection equation
- Wave equation
- Heat equation

## Discrete Integration Methods for ODEs

Before we dive into the world of PDEs, we first consider the simplest case of ODEs.

The ODE depends on a _single independent variable_ (time variable here)

We present the discrete methods for this kind of ODEs.

$$
u^{\prime}(t) = f(t, u(t))
$$

where $t \in [0, T]$, $T$ is a non-negative scalar, and $f$ is a continuous function.

**For a completely solution, we need to know the initial value of the unknonw function $u$ at $t=0$**

> This kind of problem is called _Cauchy_ or _initial value_ problem. The coupling of the ODE with an initial condition
>
> $$
> u(0) = u_0
> $$
>
> where $u_0$ is a given vector in $\mathcal{R}^m$

### Discretization

Computer can deal with a finite number of discrete values.

The numerical algorithm to solve the Cauchy problem starts by setting the points $t_0, t_1, \cdots, t_N$ at which the solution will be computed.

The points $t_n, n = 0, \cdots, N$, define a _discretization_ or a _grid_ of the interval $I = [0, T]$.

The difference between two adjacent points is called _ discretization step_, denoted with $h = t_{n+1} - t_n, n = 0, \cdots, N$.

Constant discretiztion step is the simplest case, denoted with $h = T/N$.

The sub-interval between two adjacent points is defined as $I_n = [t_n, t_{n+1}], n = 0, \cdots, N$.

> The numerical algorithm gives us a approximation of the Cauthy problem which consists in building a sequence of numbers $u_0^{(N)}, \cdots, u_N^{(N)}$ that approximate the values of $u(t_0), \cdots, u(t_N)$ of the exact solution of $u(t)$ at the same computation points.

The computation always starts with $u_0^{(N)} = u_0$ in a iterative manner, in order to satisfy the initial condition.

$u_0^{(N)}$ is a notation of vector containing the values on the grid.

## Method Based on Finite Difference

The Taylor series expansion is taken to approximate the values of the unknown $u$ for $t$ close to $t_n$.

> The discretization step $h$ being fixed, we define the following finite difference operators.

- **Forward difference**

$$
D^+u(t) = \frac{u(t+h) - u(t)}{h}
$$

- **Backward difference**

$$
D^-u(t) = \frac{u(t) - u(t - h)}{h}
$$

- **Central difference**

$$
D^0u(t) = \frac{u(t+h) - u(t-h)}{2h}
$$

### Truncation Error Denotation $\mathcal{O}$

The definition above is the result of truncation of Taylor series expansion of differentiable function $u$.

Assuming that the function $u$ is twice continuously differentiable, then there exists

$$
u(t_{n+1}) = u_t + h u^{\prime}(t_n) + \frac{h^2}{2}u^{\prime\prime}(t_n + \theta_n^+)
$$

Forward difference definition

$$
u^{\prime}(t_n) = \frac{u(t_{n=1}) - u(t_n)}{h} - \frac{h}{2}u^{\prime\prime}(t_n + \theta^+_n) \approx  D^+u(t_n)
$$

This gives us the truncation error

$$
\epsilon_n = \left| u^{\prime}(t_n) - D^+u(t_n)\right| \le \frac{h}{2}max|u^{\prime\prime}(t)|
$$

Assuming that $u^{\prime\prime}$ is bounded, then we infer that the truncation error decays to zero with $h$.

Conventionally this boundation is denoted as $\mathcal{O}(h)$

$$
u^{\prime}(t_n) = D^+u(t_n) + \mathcal{O}(h)
$$

> The order accuracy of the difference approximation is defined as the power of $h$ with which the approximation error tends to zero.

- Forward difference has first-order accuracy
- Backward difference has first-order accuracy
- Central difference has second-order accuracy

> The finite difference operators can be linearly combined to find approximations of $u^{\prime}$

$$
u^{\prime}(t_n) \approx \alpha D^-U(t_n) + \beta D^0u(t_n) + \gamma D^+ u(t_n)
$$

The coefficients $\alpha, \beta, \gamma$ are chosen such that the approximation has the highest possible order of accuray.

## Numberical Scheme for the ODE using FD approximations

- **Explicit method**

Consider the ODE at time $t_n$ and replacing $u^{\prime}(t_n)$ by $D^+ u(t_n)$, i.e., the forward difference method.

$$
u_{n+1} = u_n + hf(t_n, u_n)
$$

> The scheme is called the **explicit** Euler scheme, or simply the Euler scheme

The method is said to be explicit because $u_{n+1}$ depends explicitly on $t_n$ and the old value $u_n$

> A numerical method is said to be explicit if the unknown values can be calculated directly from quantities that are already known.

- **Implicit method**

If we replace $u^{t_{n+1}}$ by $D^-u(t_{n+1})$

$$
u_{n+1} = u_n + hf(t_{n+1}, u_{n+1})
$$

We get an _implicit_ Euler scheme this time. Computing the $u_{n+1}$ requires more work that is computing $f(t_{n+1}, u_{n+1})$

- **Leap-frog method**

If we choose the _central difference_, i.e., replacing $u^{\prime}(t_n)$ using $D^0u(t_n)$ leads to

$$
u_{n+1} = u_{n-1} + 2hf(t_n, u_n)
$$

## Numerical Integration Method

> Numerical integration is also called quadrature

Integrating the ODE on the interval $I_n = [t_n, t_{n+1}]$ yields

$$
u(t_{n+1}) - u(t_n) = \int^{t_{n+1}}_{t_n} f(s, u(s))ds = \mathcal{I}_n
$$

If we can approxmate the integral $\mathcal{I}_n$, we can compute the $u_{t_{n+1}}$ starting from the old value $u_{t_n}$.

The integral $\mathcal{I}_n$ can be apprimated through differencet ways

- **left endpoint rule**

leading to the explicit Euler scheme

$$
\mathcal{I}_n \approx hf(t_n, u_n)
$$

- **right endpoint rule**

leading to the implicit Euler scheme

$$
\mathcal{I}_n \approx hf(t_{n+1}, u_{n+1})
$$

- **midpoint rule**

$$
\mathcal{I}_n \approx hf(t_n + h/2, u(t_n + h/2))
$$

This rule leads to the modified explicit Euler scheme by using the approximation

$$
u(t_n + h/2) \approx u(t_n) + \frac{h}{2}u^{\prime}(t_n) = u(t_n) + \frac{h}{2}f(t_n, u(t_n))
$$

which gives us

$$
u_{n+1} - u_n = hf(t_n + h/2, \frac{h}{2}f(t_n, u_n))
$$

- **trapezoid rule**

leadig to the **semi-implicit** Crank-Nicolson scheme

$$
\mathcal{I}_n \approx = \frac{h}{2}\left[ f(t_n, u_n) + f(t_{n+1}, u_{n+1}) \right]
$$

## General Form of Numerical Schemes for ODEs

$$
u_{n + 1} = F(h;, t_{n+1}, u_{n+1}; t_n, u_n; \cdots)
$$

If $F$ depends on $q$ previous values $u_{u - j}, j = 0, \cdots, q-1$, the scheme is said to be a $q$-step scheme.

- The leap-frog scheme is a two-step scheme
- If $F$ does not depend on the solution at time level $t_{n+1}$, the scheme is said to be explicit, otherwise, the scheme is implicit.
