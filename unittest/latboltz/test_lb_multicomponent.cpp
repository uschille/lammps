/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "lammps.h"

#include "fix.h"
#include "info.h"
#include "modify.h"

#include "../testing/core.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstring>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {

class LbMulticomponentTest : public LAMMPSTest {
protected:
    void SetUp() override
    {
        testbinary = "LbMulticomponentTest";
        LAMMPSTest::SetUp();
    }

    // indices into the vector output of fix lb/multicomponent (inherited from lb/fluid)
    enum { TEMPERATURE = 0, MASS = 1, MOMENTUM_X = 2, MOMENTUM_Y = 3, MOMENTUM_Z = 4 };
};

// A homogeneous ternary fluid at vanishing surface tension (kappa1=kappa2=kappa3=0) is an
// exact fixed point of the lattice Boltzmann update: the free energy then reduces to the
// ideal gas expression, p0 = rho*CS2, so the equilibrium distribution of a quiescent fluid
// is feq[i] = w[i]*rho independent of the concentration fluctuations of the initial mixture.
// The update must therefore leave the total fluid mass unchanged to machine precision and
// must not generate any spurious momentum.  This pins down the mass normalization of the
// equilibrium distributions, the collide/stream step, and the halo exchange.

TEST_F(LbMulticomponentTest, HomogeneousFluidIsStationary)
{
    if (!info->has_style("fix", "lb/multicomponent")) GTEST_SKIP() << "LATBOLTZ package not built";

    BEGIN_HIDE_OUTPUT();
    command("region box block 0 16 0 16 0 8");
    command("create_box 0 box");
    command("timestep 1.0");
    command("fix lbm all lb/multicomponent 1 0.166667 1.0 D3Q19 dx 1 "
            "kappa1 0.0 kappa2 0.0 kappa3 0.0 seed 12345 init mixture");
    // the fluid moments are only available after the first update
    command("run 1 post no");
    END_HIDE_OUTPUT();

    auto *fix = lmp->modify->get_fix_by_id("lbm");
    ASSERT_NE(fix, nullptr);

    const double mass = fix->compute_vector(MASS);
    ASSERT_GT(mass, 0.0);

    BEGIN_HIDE_OUTPUT();
    command("run 10 pre no post no");
    END_HIDE_OUTPUT();

    EXPECT_DOUBLE_EQ(fix->compute_vector(MASS), mass);
    EXPECT_NEAR(fix->compute_vector(MOMENTUM_X), 0.0, 1.0e-12);
    EXPECT_NEAR(fix->compute_vector(MOMENTUM_Y), 0.0, 1.0e-12);
    EXPECT_NEAR(fix->compute_vector(MOMENTUM_Z), 0.0, 1.0e-12);
    EXPECT_NEAR(fix->compute_vector(TEMPERATURE), 0.0, 1.0e-12);
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
