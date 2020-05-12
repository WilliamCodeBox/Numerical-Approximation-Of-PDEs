#if !defined(NUMERICAL_APPXI_ODE_EULER_EXP_H)
#define NUMERICAL_APPXI_ODE_EULER_EXP_H

/**
 * @file ode-euler-exp.hpp
 * @author liu chang
 * @date 12 May 2020
 * @brief An explicit Euler scheme for ODE
 * @details
 * The methods within this file can be used to solve problems like
 * \f[
 *      u^{\prime}(t) + \alpha u(t) = f(t)
 * \f]
 *
 * The explicit Euler scheme is
 *
 * \f[
 *      u_{n+1} = u_n + hf(t_n, u_n)
 * \f]
 */

#include "types.hpp"

namespace numal {

class EulerExp {
  public:
  private:
};
} // namespace numal

#endif // NUMERICAL_APPXI_ODE_EULER_EXP_H
