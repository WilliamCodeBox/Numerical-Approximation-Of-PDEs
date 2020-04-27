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

<p align="center"><img src="/tex/6918f7123a3e8657d5634d709ff35686.svg?invert_in_darkmode&sanitize=true" align=middle width=118.63794029999998pt height=17.2895712pt/></p>

where <img src="/tex/951e2eaa1ed9be444bb82d5efe2733da.svg?invert_in_darkmode&sanitize=true" align=middle width=62.57408684999999pt height=24.65753399999998pt/>, <img src="/tex/2f118ee06d05f3c2d98361d9c30e38ce.svg?invert_in_darkmode&sanitize=true" align=middle width=11.889314249999991pt height=22.465723500000017pt/> is a non-negative scalar, and <img src="/tex/190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode&sanitize=true" align=middle width=9.81741584999999pt height=22.831056599999986pt/> is a continuous function.

**For a completely solution, we need to know the initial value of the unknonw function <img src="/tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode&sanitize=true" align=middle width=9.41027339999999pt height=14.15524440000002pt/> at <img src="/tex/1c899e1c767eb4eac89facb5d1f2cb0d.svg?invert_in_darkmode&sanitize=true" align=middle width=36.07293689999999pt height=21.18721440000001pt/>**

> This kind of problem is called _Cauchy_ or _initial value_ problem. The coupling of the ODE with an initial condition
>
> <p align="center"><img src="/tex/ba02bae3c498aabb259f7a225ea3e782.svg?invert_in_darkmode&sanitize=true" align=middle width=103.820343pt height=16.438356pt/></p>
>
> where <img src="/tex/10898c33912164da6714fe6146100886.svg?invert_in_darkmode&sanitize=true" align=middle width=15.96281939999999pt height=14.15524440000002pt/> is a given vector in <img src="/tex/45261f63bd4bf853f4e9f10c3e766610.svg?invert_in_darkmode&sanitize=true" align=middle width=25.59641699999999pt height=22.465723500000017pt/>
