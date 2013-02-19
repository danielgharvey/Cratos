#include "Particle.hpp"
#include <Exception>

template<unsigned DIM>
Particle<DIM>::Particle(unsigned index, std::vector<double> location)
	: mIndex(index),
	  mLocation(location)
{
	if (location.size() != DIM)
	{
		throw
	}
}

template<unsigned DIM>
const unsigned Particle<DIM>::GetIndex() const
{

}

template<unsigned DIM>
const std::vector<double>& Particle<DIM>::rGetLocation()
{

}
