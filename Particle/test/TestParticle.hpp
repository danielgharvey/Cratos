// Test Particle

#include <cxxtest/TestSuite.h>
#include "Particle.hpp"

class TestParticle : public CxxTest::TestSuite
{
public:

    void TestConstructor()
    {
    	std::vector<double> location(3,0);

    	TS_ASSERT_THROWS(Particle<1> particle(0, location), std::runtime_error);
    	TS_ASSERT_THROWS(Particle<1> particle(0, location), std::runtime_error);

    	Particle<3> particle(0, location);

    	TS_ASSERT_EQUALS(particle.rGetLocation().size(), 3u);
    	TS_ASSERT_DELTA(particle.rGetLocation()[0], 0.0, 1e-10);
    	TS_ASSERT_DELTA(particle.rGetLocation()[1], 0.0, 1e-10);
    	TS_ASSERT_DELTA(particle.rGetLocation()[2], 0.0, 1e-10);

    	TS_ASSERT_EQUALS(particle.GetIndex(), 0u);
    }
};
