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
