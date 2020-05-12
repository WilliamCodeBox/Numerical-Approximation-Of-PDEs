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

### Discretization

Computer can deal with a finite number of discrete values.

The numerical algorithm to solve the Cauchy problem starts by setting the points <img src="/tex/2068048116354b91f18e504c12d02907.svg?invert_in_darkmode&sanitize=true" align=middle width=90.77833215pt height=20.221802699999984pt/> at which the solution will be computed.

The points <img src="/tex/91069424bdc68b03544149ad16c40901.svg?invert_in_darkmode&sanitize=true" align=middle width=116.46270075pt height=22.465723500000017pt/>, define a _discretization_ or a _grid_ of the interval <img src="/tex/6d452ef266247489a5a92bfe866693ad.svg?invert_in_darkmode&sanitize=true" align=middle width=66.98044814999999pt height=24.65753399999998pt/>.

The difference between two adjacent points is called _ discretization step_, denoted with <img src="/tex/2fa709a95cab1683601c145e9f2f45aa.svg?invert_in_darkmode&sanitize=true" align=middle width=199.4705823pt height=22.831056599999986pt/>.

Constant discretiztion step is the simplest case, denoted with <img src="/tex/ff0c3443727e003a74efe7ad4c273a07.svg?invert_in_darkmode&sanitize=true" align=middle width=65.12737604999998pt height=24.65753399999998pt/>.

The sub-interval between two adjacent points is defined as <img src="/tex/3a9ac8cfc4e36fea359fedf5a64284e5.svg?invert_in_darkmode&sanitize=true" align=middle width=202.52059574999996pt height=24.65753399999998pt/>.

> The numerical algorithm gives us a approximation of the Cauthy problem which consists in building a sequence of numbers <img src="/tex/e91099538bdbdc9c58fd266677e1a66b.svg?invert_in_darkmode&sanitize=true" align=middle width=102.75186734999997pt height=34.337843099999986pt/> that approximate the values of <img src="/tex/6edc1ffd5c683540010e975d677d01e9.svg?invert_in_darkmode&sanitize=true" align=middle width=115.37519564999998pt height=24.65753399999998pt/> of the exact solution of <img src="/tex/633aafb63e9ef3f6da5d763d63ed3a95.svg?invert_in_darkmode&sanitize=true" align=middle width=28.13180369999999pt height=24.65753399999998pt/> at the same computation points.

The computation always starts with <img src="/tex/231c350084355c561e9a3197afc18d29.svg?invert_in_darkmode&sanitize=true" align=middle width=70.03279799999999pt height=34.337843099999986pt/> in a iterative manner, in order to satisfy the initial condition.

<img src="/tex/daf8d091e66a28ae10625dd2dcb0e53e.svg?invert_in_darkmode&sanitize=true" align=middle width=31.33044419999999pt height=34.337843099999986pt/> is a notation of vector containing the values on the grid.

## Method Based on Finite Difference

The Taylor series expansion is taken to approximate the values of the unknown <img src="/tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode&sanitize=true" align=middle width=9.41027339999999pt height=14.15524440000002pt/> for <img src="/tex/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode&sanitize=true" align=middle width=5.936097749999991pt height=20.221802699999984pt/> close to <img src="/tex/27413cd33c6f718117d8fb364284f787.svg?invert_in_darkmode&sanitize=true" align=middle width=14.06212004999999pt height=20.221802699999984pt/>.

> The discretization step <img src="/tex/2ad9d098b937e46f9f58968551adac57.svg?invert_in_darkmode&sanitize=true" align=middle width=9.47111549999999pt height=22.831056599999986pt/> being fixed, we define the following finite difference operators.

- **Forward difference**

<p align="center"><img src="/tex/d955713beedea9228e0f3be3b2ccf12a.svg?invert_in_darkmode&sanitize=true" align=middle width=182.9186535pt height=34.7253258pt/></p>

- **Backward difference**

<p align="center"><img src="/tex/1d47905c8dbdd49e9752185b7a6a6d4f.svg?invert_in_darkmode&sanitize=true" align=middle width=183.10128375pt height=34.7253258pt/></p>

- **Central difference**

<p align="center"><img src="/tex/15cc18d19b65aa15a848380ef9115a0f.svg?invert_in_darkmode&sanitize=true" align=middle width=208.9421334pt height=34.7253258pt/></p>

### Truncation Error Denotation <img src="/tex/9fa4bf66c871f8af69c9d3cf2fcb6a55.svg?invert_in_darkmode&sanitize=true" align=middle width=13.54343924999999pt height=22.465723500000017pt/>

The definition above is the result of truncation of Taylor series expansion of differentiable function <img src="/tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode&sanitize=true" align=middle width=9.41027339999999pt height=14.15524440000002pt/>.

Assuming that the function <img src="/tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode&sanitize=true" align=middle width=9.41027339999999pt height=14.15524440000002pt/> is twice continuously differentiable, then there exists

<p align="center"><img src="/tex/76f95cd3cdf82b36198b08f453544c1e.svg?invert_in_darkmode&sanitize=true" align=middle width=287.63468415pt height=35.77743345pt/></p>

Forward difference definition

<p align="center"><img src="/tex/c73fc67029ca32b1f777b4384f1799c7.svg?invert_in_darkmode&sanitize=true" align=middle width=380.59290719999996pt height=34.7253258pt/></p>

This gives us the truncation error

<p align="center"><img src="/tex/96be73aff2aa89178fa41c5079401644.svg?invert_in_darkmode&sanitize=true" align=middle width=285.8560716pt height=33.81208709999999pt/></p>

Assuming that <img src="/tex/3ca6d260772e6d85c100c7e3f4888e48.svg?invert_in_darkmode&sanitize=true" align=middle width=16.99019684999999pt height=24.7161288pt/> is bounded, then we infer that the truncation error decays to zero with <img src="/tex/2ad9d098b937e46f9f58968551adac57.svg?invert_in_darkmode&sanitize=true" align=middle width=9.47111549999999pt height=22.831056599999986pt/>.

Conventionally this boundation is denoted as <img src="/tex/34109b622bca0e0d13a8a0d0cca985dd.svg?invert_in_darkmode&sanitize=true" align=middle width=35.79998729999999pt height=24.65753399999998pt/>

<p align="center"><img src="/tex/80945c83248da22578b0d551637ccd04.svg?invert_in_darkmode&sanitize=true" align=middle width=181.55968545pt height=18.0201615pt/></p>

> The order accuracy of the difference approximation is defined as the power of <img src="/tex/2ad9d098b937e46f9f58968551adac57.svg?invert_in_darkmode&sanitize=true" align=middle width=9.47111549999999pt height=22.831056599999986pt/> with which the approximation error tends to zero.

- Forward difference has first-order accuracy
- Backward difference has first-order accuracy
- Central difference has second-order accuracy

> The finite difference operators can be linearly combined to find approximations of <img src="/tex/ac2b68ee3af718535820a82765ee436f.svg?invert_in_darkmode&sanitize=true" align=middle width=13.200234299999991pt height=24.7161288pt/>

<p align="center"><img src="/tex/73618d796afe3c3c299960171691b641.svg?invert_in_darkmode&sanitize=true" align=middle width=320.38476855pt height=18.312383099999998pt/></p>

The coefficients <img src="/tex/d3bd6dd6baaf64c04cdf54776ebf0b37.svg?invert_in_darkmode&sanitize=true" align=middle width=44.777676899999996pt height=22.831056599999986pt/> are chosen such that the approximation has the highest possible order of accuray.

## Numberical Scheme for the ODE using FD approximations

- **Explicit method**

Consider the ODE at time <img src="/tex/27413cd33c6f718117d8fb364284f787.svg?invert_in_darkmode&sanitize=true" align=middle width=14.06212004999999pt height=20.221802699999984pt/> and replacing <img src="/tex/057215554dc669e2f4df335b36d425c5.svg?invert_in_darkmode&sanitize=true" align=middle width=41.69161589999999pt height=24.7161288pt/> by <img src="/tex/701797be0c0c29838f70d95a54a69d49.svg?invert_in_darkmode&sanitize=true" align=middle width=62.05925879999998pt height=26.17730939999998pt/>, i.e., the forward difference method.

<p align="center"><img src="/tex/861d0e5140d1259e7f957ed9b61f5bf5.svg?invert_in_darkmode&sanitize=true" align=middle width=167.99124375pt height=16.438356pt/></p>

> The scheme is called the **explicit** Euler scheme, or simply the Euler scheme

The method is said to be explicit because <img src="/tex/2a05f4d7c156553a40d8fc307a37ed49.svg?invert_in_darkmode&sanitize=true" align=middle width=34.180216949999995pt height=14.15524440000002pt/> depends explicitly on <img src="/tex/27413cd33c6f718117d8fb364284f787.svg?invert_in_darkmode&sanitize=true" align=middle width=14.06212004999999pt height=20.221802699999984pt/> and the old value <img src="/tex/f38fb14fa5981bc6b0514cebd4cd10ab.svg?invert_in_darkmode&sanitize=true" align=middle width=17.53629569999999pt height=14.15524440000002pt/>

> A numerical method is said to be explicit if the unknown values can be calculated directly from quantities that are already known.

- **Implicit method**

If we replace <img src="/tex/6980ef95a908041256d1c7ab676d3ee5.svg?invert_in_darkmode&sanitize=true" align=middle width=35.65686299999999pt height=26.085962100000025pt/> by <img src="/tex/cced5776bd4b22bbcd237e65b4799ab4.svg?invert_in_darkmode&sanitize=true" align=middle width=78.88579874999999pt height=26.17730939999998pt/>

<p align="center"><img src="/tex/eff46a79834c33762b392f4c840a90a4.svg?invert_in_darkmode&sanitize=true" align=middle width=201.27906645pt height=16.438356pt/></p>

We get an _implicit_ Euler scheme this time. Computing the <img src="/tex/2a05f4d7c156553a40d8fc307a37ed49.svg?invert_in_darkmode&sanitize=true" align=middle width=34.180216949999995pt height=14.15524440000002pt/> requires more work that is computing <img src="/tex/4a0dfdfa57473cbf3ad81c19fc9c8bcc.svg?invert_in_darkmode&sanitize=true" align=middle width=96.43879739999998pt height=24.65753399999998pt/>

- **Leap-frog method**

If we choose the _central difference_, i.e., replacing <img src="/tex/057215554dc669e2f4df335b36d425c5.svg?invert_in_darkmode&sanitize=true" align=middle width=41.69161589999999pt height=24.7161288pt/> using <img src="/tex/b217407422390158e850ea2b4b608001.svg?invert_in_darkmode&sanitize=true" align=middle width=58.520432849999985pt height=26.76175259999998pt/> leads to

<p align="center"><img src="/tex/518c219ab9ff594d54865edddfc80220.svg?invert_in_darkmode&sanitize=true" align=middle width=193.03699304999998pt height=16.438356pt/></p>

## Numerical Integration Method

> Numerical integration is also called quadrature

Integrating the ODE on the interval <img src="/tex/64364780ee67c723b87f2282b1566482.svg?invert_in_darkmode&sanitize=true" align=middle width=100.94193119999997pt height=24.65753399999998pt/> yields

<p align="center"><img src="/tex/a1ef8171487f2e4f17d3195217e4245a.svg?invert_in_darkmode&sanitize=true" align=middle width=302.65206344999996pt height=42.00953955pt/></p>

If we can approxmate the integral <img src="/tex/779c7b09279e8270862dc891ff4c85a8.svg?invert_in_darkmode&sanitize=true" align=middle width=17.07696539999999pt height=22.465723500000017pt/>, we can compute the <img src="/tex/358d269de47a5afe6f87e3a3f587ebb7.svg?invert_in_darkmode&sanitize=true" align=middle width=35.65686299999999pt height=14.15524440000002pt/> starting from the old value <img src="/tex/a895d08ccd87e48eb3754d811de09dfa.svg?invert_in_darkmode&sanitize=true" align=middle width=21.615526349999993pt height=14.15524440000002pt/>.

The integral <img src="/tex/779c7b09279e8270862dc891ff4c85a8.svg?invert_in_darkmode&sanitize=true" align=middle width=17.07696539999999pt height=22.465723500000017pt/> can be apprimated through differencet ways

- **left endpoint rule**

leading to the explicit Euler scheme

<p align="center"><img src="/tex/07e9c9fbc3c406b690108af01f06fde9.svg?invert_in_darkmode&sanitize=true" align=middle width=112.43859989999999pt height=16.438356pt/></p>

- **right endpoint rule**

leading to the implicit Euler scheme

<p align="center"><img src="/tex/d1e4871baedfc9af10ae7551b4d21756.svg?invert_in_darkmode&sanitize=true" align=middle width=145.7264226pt height=16.438356pt/></p>

- **midpoint rule**

<p align="center"><img src="/tex/9b1fce5210ef2920b9789ab14423ed36.svg?invert_in_darkmode&sanitize=true" align=middle width=223.16158095pt height=16.438356pt/></p>

This rule leads to the modified explicit Euler scheme by using the approximation

<p align="center"><img src="/tex/52fae16975ccb6f5e9f0ed433c6c0ac7.svg?invert_in_darkmode&sanitize=true" align=middle width=391.6542927pt height=33.81208709999999pt/></p>

which gives us

<p align="center"><img src="/tex/d49fb7c80a5fdc42599bac7d142ca5c1.svg?invert_in_darkmode&sanitize=true" align=middle width=272.2010214pt height=33.81208709999999pt/></p>

- **trapezoid rule**

leadig to the **semi-implicit** Crank-Nicolson scheme

<p align="center"><img src="/tex/bd41012b1ff1c2de4211e33e10321967.svg?invert_in_darkmode&sanitize=true" align=middle width=257.571303pt height=33.81208709999999pt/></p>

## General Form of Numerical Schemes for ODEs

<p align="center"><img src="/tex/ae502f765ddde50d94239ff29e9adf11.svg?invert_in_darkmode&sanitize=true" align=middle width=257.55548445pt height=16.438356pt/></p>

If <img src="/tex/b8bc815b5e9d5177af01fd4d3d3c2f10.svg?invert_in_darkmode&sanitize=true" align=middle width=12.85392569999999pt height=22.465723500000017pt/> depends on <img src="/tex/d5c18a8ca1894fd3a7d25f242cbe8890.svg?invert_in_darkmode&sanitize=true" align=middle width=7.928106449999989pt height=14.15524440000002pt/> previous values <img src="/tex/77d35a1cc6ebb1f8cb2be596715ca6e2.svg?invert_in_darkmode&sanitize=true" align=middle width=155.04352709999998pt height=21.68300969999999pt/>, the scheme is said to be a <img src="/tex/d5c18a8ca1894fd3a7d25f242cbe8890.svg?invert_in_darkmode&sanitize=true" align=middle width=7.928106449999989pt height=14.15524440000002pt/>-step scheme.

- The leap-frog scheme is a two-step scheme
- If <img src="/tex/b8bc815b5e9d5177af01fd4d3d3c2f10.svg?invert_in_darkmode&sanitize=true" align=middle width=12.85392569999999pt height=22.465723500000017pt/> does not depend on the solution at time level <img src="/tex/a9b85307359e3224805c2d3b5192873a.svg?invert_in_darkmode&sanitize=true" align=middle width=30.70604129999999pt height=20.221802699999984pt/>, the scheme is said to be explicit, otherwise, the scheme is implicit.
