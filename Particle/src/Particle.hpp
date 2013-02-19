/*
 * A general particle class.
 */
template<unsigned DIM>
class Particle
{
private:

  unsigned mIndex;

public:

  Particle()
  {
    mIndex = 0;
  }

  unsigned GetIndex()
  {
    return mIndex;
  }
};
