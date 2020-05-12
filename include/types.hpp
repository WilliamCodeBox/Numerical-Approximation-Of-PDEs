#if !defined(NUMALALGOPDES_TYPES_H)
#define NUMALALGOPDES_TYPES_H

#include <Eigen/Core>
#include <complex>

namespace numal {
/*short hand for std double complex type */
using dcomplex = std::complex<double>;

/* short hand for double type Eigen matrix sized 3x1 */
using vec3 = Eigen::Matrix<double, 3, 1>;

/* short hand for double complex type Eigen matrix sized 3x1 */
using cvec3 = Eigen::Matrix<dcomplex, 3, 1>;
} // namespace numal

#endif // NUMALALGOPDES_TYPES_H
