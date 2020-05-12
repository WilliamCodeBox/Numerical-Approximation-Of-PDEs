#include "types.hpp"
#include <catch2/catch.hpp>

namespace numal {
TEST_CASE("Basic data types", "[numal, Types]") {

    SECTION("double complex type") {
        dcomplex a{1, 2};
        REQUIRE(a.real() == 1.0);
        REQUIRE(a.imag() == 2.0);
    }
}
} // namespace numal
