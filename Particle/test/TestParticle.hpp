// Test Particle

#include <cxxtest/TestSuite.h>
#include "Particle.hpp"

class TestParticle : public CxxTest::TestSuite
{
public:

    void TestConstructor()
    {
        Particle<1> particle;
        TS_ASSERT_EQUALS(particle.GetIndex(), 1u);
    }
};
