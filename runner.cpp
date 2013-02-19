/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#define _CXXTEST_HAVE_EH
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/TestMain.h>
#include <cxxtest/ErrorPrinter.h>

int main( int argc, char *argv[] ) {
 int status;
    CxxTest::ErrorPrinter tmp;
    status = CxxTest::Main<CxxTest::ErrorPrinter>( tmp, argc, argv );
    return status;
}
bool suite_TestParticle_init = false;
#include "Particle/test/TestParticle.hpp"

static TestParticle suite_TestParticle;

static CxxTest::List Tests_TestParticle = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestParticle( "Particle/test/TestParticle.hpp", 6, "TestParticle", suite_TestParticle, Tests_TestParticle );

static class TestDescription_suite_TestParticle_TestConstructor : public CxxTest::RealTestDescription {
public:
 TestDescription_suite_TestParticle_TestConstructor() : CxxTest::RealTestDescription( Tests_TestParticle, suiteDescription_TestParticle, 10, "TestConstructor" ) {}
 void runTest() { suite_TestParticle.TestConstructor(); }
} testDescription_suite_TestParticle_TestConstructor;

#include <cxxtest/Root.cpp>
const char* CxxTest::RealWorldDescription::_worldName = "cxxtest";
