/*
 * A general particle class.
 */
template<unsigned DIM>
class Particle
{
private:

	std::vector<double> mLocation;

	unsigned mIndex;

public:

	Particle(unsigned index, std::vector<double> location);

	const unsigned GetIndex() const;

	const std::vector<double>& rGetLocation();
};
