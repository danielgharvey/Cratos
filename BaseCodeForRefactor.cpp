/*
* Parallel implementation of simple particle simulator for particles on a square domain with reflective boundary conditions. Neighbouring particles are assumed to
* be connected by a linear spring. Random motion is implemented using a 'random force'.
*/

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <mpi.h>
#include <time.h>
#include <set>
#include <assert.h>
#include <map>

using namespace std;


/* A simplified particle */
class Particle
{
private:
    /* Location of the particle centre */
    vector<double> mLocation;

    /* Index of the particle */
    unsigned mIndex;

public:

    /**
     * Constructor
     *
     * @param rLocation the locaiton of the particle centre.
     */
    Particle(unsigned index, vector<double>& rLocation)
    {
	assert(rLocation.size() ==  2);
	mIndex = index;
	mLocation = rLocation;
    }

    unsigned GetIndex()
    {
	return mIndex;
    }

    std::vector<double>& rGetLocation()
    {
	return mLocation;
    }

};

/* A simplified version of the box class in Chaste */
class Box
{
private:
    /** Coordinates of the box, in the form (for 2D) (xmin, xmax, ymin, ymax) (etc). */
    std::vector<double> mMinAndMaxValues;

    /** Nodes contained in this box. */
    std::set<Particle*> mParticlesContained;

public:

    /**
     * Constructor just takes in the extremal values of the box.
     *
     * @param rMinAndMaxValues the extremal values. Of the from (for 2D, etc): xmin, xmax, ymin, ymax
     */
    Box(vector<double>& rMinAndMaxValues)
    {
	    mMinAndMaxValues = rMinAndMaxValues;
    }

    /** Get the coordinates of the box, in the form (for 2D) (xmin, xmax, ymin, ymax) (etc). */
    vector<double>& rGetMinAndMaxValues()
    {
	return mMinAndMaxValues;
    }

    /**
     * Add a node to this box.
     * @param pNode address of the node to be added
     */
    void AddParticle(Particle* p_Particle)
    {
	mParticlesContained.insert(p_Particle);
    }

    /**
     * Remove a node from this box.
     * @param pNode address of the node to be removed
     */
    void RemoveParticle(Particle* p_Particle)
    {
    	mParticlesContained.erase(p_Particle);
    }

    /**
     * Clear all the Particles from the box
     */
    void Clear()
    {
	mParticlesContained.clear();
    }

    /** Get all the nodes in this box. */
    std::set<Particle*>& rGetParticlesContained()
    {
	return mParticlesContained;
    }

    /** Is point contained in box */
    bool IsPointContained(std::vector<double>& rPoint)
    {
	assert(rPoint.size() == 2);
	std::cout << "(" << rPoint[0] << "," << rPoint[1] << ")" << "\n";
	return ((rPoint[0]>mMinAndMaxValues[0] || rPoint[0] == mMinAndMaxValues[0]) && (rPoint[0] < mMinAndMaxValues[1]) && (rPoint[1] > mMinAndMaxValues[2] || rPoint[1] ==  mMinAndMaxValues[2]) && (rPoint[1] < mMinAndMaxValues[3]));
    }

};

unsigned VectorMin(std::vector<double> vec)
{
	//Assume it is zero
	unsigned min = 0;

	//Test otherwise
	for (unsigned i = 0; i<vec.size(); i++)
	{
		min = (vec[i]<vec[min]) ? i : min;
	}

	return min;
}

unsigned MinimumMetric(double x, double y, double domain_size)
{
	double xfactor = (x>0) ? domain_size : (-1.0)*domain_size;	
	double yfactor = (y>0) ? domain_size : (-1.0)*domain_size;

	// Different metrics to account for periodicity
	std::vector<double> d(4,0.0);
	d[0] = x*x+y*y;
	d[1] = (x-xfactor)*(x-xfactor)+y*y;
	d[2] = x*x+(y-yfactor)*(y-yfactor);
	d[3] = (x-xfactor)*(x-xfactor)+(y-yfactor)*(y-yfactor);

	// Find the minimum
	unsigned min_index = VectorMin(d);

	return min_index;
}

double DistanceUsingMetric(double x, double y, double domain_size, unsigned metric)
{
	double xfactor = (x>0) ? domain_size : (-1.0)*domain_size;	
	double yfactor = (y>0) ? domain_size : (-1.0)*domain_size;

	// Different metrics to account for periodicity
	std::vector<double> d(4,0.0);
	d[0] = x*x+y*y;
	d[1] = (x-xfactor)*(x-xfactor)+y*y;
	d[2] = x*x+(y-yfactor)*(y-yfactor);
	d[3] = (x-xfactor)*(x-xfactor)+(y-yfactor)*(y-yfactor);

	// Find the minimum
	unsigned min_index = VectorMin(d);

	return d[metric];
}


/*
* Given a set of particles and a distance, determine pairs of particles that lie within that distance
*/
std::vector<std::pair<unsigned, unsigned> > FindPairs(std::vector< std::vector<double> > particles, double mechanicslength, double domain_size)
{
	// Allocate memory
	std::vector<std::pair<unsigned, unsigned> > Pairs;
	
	for (unsigned i = 0; i < particles.size(); i++)
	{
		for (unsigned j = 0; j < particles.size(); j++)
		{

			// Calculate distance including periodicity conditions.
			double x = particles[i][0]-particles[j][0];
			double y = particles[i][1]-particles[j][1];
	
			double min_index = MinimumMetric(x,y,domain_size);
			double d = DistanceUsingMetric(x,y,domain_size,min_index);
				
			// Add close particles to list of pairs. First clause here prevent duplicate pairs being recorded
			if ( d > 0.000001 && d < mechanicslength)
			{
				std::pair<unsigned, unsigned> temp;
				temp.first = i;
				temp.second = j;
				Pairs.push_back(temp);
			}
		}
	}
	return Pairs;		
}


std::vector<std::vector<double> > CalculateRepulsionForce(std::vector< std::vector<double> > particles, double mechanicslength, double domain_size)
{
	// Find Pairs

	std::vector<std::pair<unsigned, unsigned> > Pairs;
	
	Pairs = FindPairs(particles, mechanicslength, domain_size);

	// Create memory to store foce contributions. Initialised to zero for isolated particles.
	std::vector<std::vector<double> > ForceContribution(particles.size(), std::vector<double>(2,0));	

	for (unsigned i = 0; i < Pairs.size(); i++)
	{
		unsigned NodeAIndex = Pairs[i].first;
		unsigned NodeBIndex = Pairs[i].second;

		std::vector<double> Separation(2,0);
	
		// Calculate separating vector
		Separation[0] = particles[NodeBIndex][0]-particles[NodeAIndex][0];
		Separation[1] = particles[NodeBIndex][1]-particles[NodeAIndex][1];

		// Establish correct metric

		std::vector<double> d(4,0);
		double x = Separation[0];
		double y = Separation[1];

		unsigned metric = MinimumMetric(x,y,domain_size);
		double Distance = DistanceUsingMetric(x,y,domain_size,metric);

		// Get unit vector connecting particles (we take into account periodic direction in calculating the force).

		double unit = Separation[0]*Separation[0]+Separation[1]*Separation[1];
		Separation[0] = Separation[0]/unit;
		Separation[1] = Separation[1]/unit;

		// Calculate magnitude of force between particles and add force contribution to vector of forces. See Meineke 2001 (DOI: 10.1046/j.0960-7722.2001.00216.x) for details.

		double SpringConstant = 0.1;
		double SpringRestLength = 0.9;

		switch(metric)
		{
			case 0:

				ForceContribution[NodeAIndex][0]+= (-1)*SpringConstant*(SpringRestLength-Distance)*Separation[0];
				ForceContribution[NodeAIndex][1]+= (-1)*SpringConstant*(SpringRestLength-Distance)*Separation[1];
				ForceContribution[NodeBIndex][0]+= SpringConstant*(SpringRestLength-Distance)*Separation[0];
				ForceContribution[NodeBIndex][1]+= SpringConstant*(SpringRestLength-Distance)*Separation[1];
				
				break;

			case 1:
				ForceContribution[NodeAIndex][0]+= SpringConstant*(SpringRestLength-Distance)*Separation[0];
				ForceContribution[NodeAIndex][1]+= (-1)*SpringConstant*(SpringRestLength-Distance)*Separation[1];
				ForceContribution[NodeBIndex][0]+= (-1)*SpringConstant*(SpringRestLength-Distance)*Separation[0];
				ForceContribution[NodeBIndex][1]+= SpringConstant*(SpringRestLength-Distance)*Separation[1];
				
				break;		
			case 2:
				ForceContribution[NodeAIndex][0]+= (-1)*SpringConstant*(SpringRestLength-Distance)*Separation[0];
				ForceContribution[NodeAIndex][1]+= SpringConstant*(SpringRestLength-Distance)*Separation[1];
				ForceContribution[NodeBIndex][0]+= SpringConstant*(SpringRestLength-Distance)*Separation[0];
				ForceContribution[NodeBIndex][1]+= (-1)*SpringConstant*(SpringRestLength-Distance)*Separation[1];
				
				break;
			case 3:
				ForceContribution[NodeAIndex][0]+= SpringConstant*(SpringRestLength-Distance)*Separation[0];
				ForceContribution[NodeAIndex][1]+= SpringConstant*(SpringRestLength-Distance)*Separation[1];
				ForceContribution[NodeBIndex][0]+= (-1)*SpringConstant*(SpringRestLength-Distance)*Separation[0];
				ForceContribution[NodeBIndex][1]+= (-1)*SpringConstant*(SpringRestLength-Distance)*Separation[1];
				
				break;
		}	
	}

	return ForceContribution;
}

std::vector< std::vector<double> > CalculateRandomForce(std::vector< Particle* >& rParticles, double nu, double dt, double D)
{
	// Create force contribution vector
	std::vector<std::vector<double> > ForceContribution(rParticles.size(), std::vector<double>(2,0));

	for (unsigned i = 0; i < rParticles.size(); i++)
	{

		// Generate random direction
		double u = (double)rand()/RAND_MAX - 1.0;
		double v = (double)rand()/RAND_MAX - 1.0;

		double mag = sqrt(u*u + v*v);

		double X = u / mag;
		double Y = v / mag;

		// Add force contribution
		ForceContribution[i][0] = (nu*sqrt(2.0*D*dt)/dt)*X;
		ForceContribution[i][1] = (nu*sqrt(2.0*D*dt)/dt)*Y;

	}
	return ForceContribution;
}

/*
 * Calculate the global index of the box that contains a particle at a given point
 * @param ParticleCentre the location of the particle centre.
 */

unsigned CalculateContainingBox(std::vector<double> ParticleCentre)
{
	assert(ParticleCentre.size() ==  2);
	// THIS MEHTOD MAKES ASSUMPTION THAT DOMAIN SIZE IS 300 AND BOX WIDTH IS 1.5. IF THIS CHANGES THIS METHOD SHOULD BE CHANGED TO REFLECT THIS	
	unsigned location;
	
	unsigned x_comp = floor (ParticleCentre[0] / 1.5);
	unsigned y_comp = floor (ParticleCentre[1] / 1.5);

	location = x_comp + 200*y_comp;

	return location;	
}

/*
 * Return the component of the vector (0,1,2,...199) owned by this process. For use in setting up the domain
 * @param pid the id of this process
 * @param numprocs the total number of processes.
 */
std::vector<unsigned> GetDistributedVectorComponent(unsigned pid, unsigned numprocs)
{
	unsigned vec_size;
	if (pid < numprocs-1)
	{
		vec_size = ceil(200.0/numprocs);
	}
	else
	{
		vec_size = 200 - (numprocs-1) * ceil (200.0/numprocs);
	}
	
	std::vector<unsigned> vec(vec_size,0);
	
	// Populate vector
	for (unsigned i = 0; i<vec.size(); i++)
	{
		vec[i] = pid * ceil(200.0/numprocs) + i;
	}

	return vec;
}

bool IsParticleOwned(Particle* pParticle, std::map<unsigned, unsigned>& rBoxesGlobalToLocalIndexMapping)
{
	bool is_owned;
	
	// Calculate the containing box
	unsigned containing_box = CalculateContainingBox(pParticle->rGetLocation());
	
	return (rBoxesGlobalToLocalIndexMapping.find(containing_box) !=  rBoxesGlobalToLocalIndexMapping.end());	
}


int main(int argv, char **argc)
{
	// Initialise MPI
	int numprocs, pid, right, left;
	MPI_Init(&argv, &argc);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);	
	MPI_Comm_rank(MPI_COMM_WORLD,&pid);
	MPI_Status status;
	
	// Whether to write output to file. Set to 'false' for benchmarking tests.
	bool WriteOutputToFiles = false;
	
	// Name file according to NP for easy analysis
	stringstream ss (stringstream::in | stringstream::out);
	string np, filename;
	ss << numprocs;

	np = ss.str();
	filename = "OutputNP"+np+".txt";

	// Open file for results
	std::ofstream file(filename.c_str());


	// Timing Variables
	double begin, end, mpi_time;
	double start = MPI_Wtime();
	double time_calculating_pairs = 0.0; 

	right =  (pid + 1) % numprocs;
    	left = pid - 1;
    	if (left < 0)
	{
        	left = numprocs - 1;
	}

	bool AmMaster = (pid ==  0);
	bool AmTopMost = (pid ==  numprocs -1);
	
	/*
	* Set up local memory
	*/

	//Parameters
	unsigned TimeSteps = 120;
	double EndTime = 1.0;
	double DeltaT = EndTime/(TimeSteps-1);
	double DomainSize = 300.0;
	double DampingConstant = 0.1; 
	double MechanicsCutOff = 1.5;

	/*************************
	 * Set up the boxes
	************************/
	std::vector<unsigned> vector_comp = GetDistributedVectorComponent(pid, numprocs);


	// Set up bound of the process [xmin, xmax)
	std::vector<double> RegionBound(2,0);
	RegionBound[0] = (double)(1.5 * vector_comp[0]);
	RegionBound[1] = (double)(1.5 * (vector_comp[vector_comp.size() - 1] + 1));

	// Check that we aren't using too many processes for the domain size
	if (RegionBound[1] - RegionBound[0] < MechanicsCutOff)
	{
		std::cout << "Too many processes in use for this domain size and mechanics cut off.\n";
	}

	// Create boxes
	std::vector<Box> mBoxes;

	// A map from global to local box indices
	std::map<unsigned, unsigned> mBoxGlobalToLocalIndexMapping;
	// And the inverse map.
	std::map<unsigned, unsigned> mBoxLocalToGlobalIndexMapping;

	for (unsigned y_comp = 0; y_comp<200; y_comp++)
	{	
		for (unsigned i = 0; i<vector_comp.size(); i++)
		{
			std::vector<double> box_location(4,0);
			box_location[0] = 1.5 * vector_comp[i];
			box_location[1] = 1.5 * (vector_comp[i]+1);
			box_location[2] = 1.5 * y_comp;
			box_location[3] = 1.5 * (y_comp + 1);

			Box new_box(box_location);
			mBoxes.push_back(new_box);

			// Update the map between global and local box indices.
			unsigned global_index = vector_comp[i] + 200 * y_comp;
			mBoxGlobalToLocalIndexMapping[global_index] = mBoxes.size()-1;
			mBoxLocalToGlobalIndexMapping[mBoxes.size()-1] = global_index;
		}
	}

	// Set up the local boxes topology - based on global indices.
	std::vector< std::set<unsigned> > mLocalBoxes;

	unsigned M = 200;
	unsigned N = 200;

	std::vector<bool> is_xmin(N*M); // far left
	std::vector<bool> is_xmax(N*M); // far right
	std::vector<bool> is_ymin(N*M); // bottom
	std::vector<bool> is_ymax(N*M); // top

	for (unsigned i = 0; i<M*N; i++)
	{
		is_xmin[i] = (i%M == 0);
		is_xmax[i] = ((i+1)%M == 0);
		is_ymin[i] = (i%(M*N)<M);
		is_ymax[i] = (i%(M*N)>= (N-1)*M);
	}

	for (unsigned local_index = 0; local_index < mBoxes.size(); local_index++)
	{
		// Set i to be the global index
		unsigned i = mBoxLocalToGlobalIndexMapping[local_index];

		std::set<unsigned> local_boxes;

		local_boxes.insert(i);

		// add the box to the left
		if (!is_xmin[i])
		{
			local_boxes.insert(i-1);
		}

		// add the box to the right
		if (!is_xmax[i])
		{
			local_boxes.insert(i+1);
		}

		// add the one below
		if (!is_ymin[i])
		{
			local_boxes.insert(i-M);
		}

		// add the one above
		if (!is_ymax[i])
		{
			local_boxes.insert(i+M);
		}

		// add the four corner boxes

		if ( (!is_xmin[i]) && (!is_ymin[i]) )
		{
			local_boxes.insert(i-1-M);
		}

		if ( (!is_xmin[i]) && (!is_ymax[i]) )
		{
			local_boxes.insert(i-1+M);
		}

		if ( (!is_xmax[i]) && (!is_ymin[i]) )
		{
			local_boxes.insert(i+1-M);
		}

		if ( (!is_xmax[i]) && (!is_ymax[i]) )
		{
			local_boxes.insert(i+1+M);
		}

		mLocalBoxes.push_back(local_boxes);
	}

	// Create halo boxes on each process.
	std::vector<Box> mHaloBoxesRight;
	std::vector<Box> mHaloBoxesLeft;

	// A map from global to local halo box indices
	std::map<unsigned, unsigned> mHaloBoxRightGlobalToLocalIndexMapping;
	std::map<unsigned, unsigned> mHaloBoxLeftGlobalToLocalIndexMapping;
	// And the inverse map.
	std::map<unsigned, unsigned> mHaloBoxRightLocalToGlobalIndexMapping;
	std::map<unsigned, unsigned> mHaloBoxLeftLocalToGlobalIndexMapping;

	for (unsigned y_comp = 0; y_comp<200; y_comp++)
	{
		if (!AmTopMost)
		{
			std::vector<double> box_location(4,0);
			box_location[0] = 1.5 * (vector_comp[vector_comp.size()-1]+1);
			box_location[1] = 1.5 * (vector_comp[vector_comp.size()-1]+2);
			box_location[2] = 1.5 * y_comp;
			box_location[3] = 1.5 * (y_comp + 1);

			Box new_box(box_location);
			mHaloBoxesRight.push_back(new_box);

			// Update the map between global and local box indices.
			unsigned global_index = vector_comp[vector_comp.size()-1] + 1 + 200 * y_comp;
			mHaloBoxRightGlobalToLocalIndexMapping[global_index] = mHaloBoxesRight.size()-1;
			mHaloBoxRightLocalToGlobalIndexMapping[mHaloBoxesRight.size()-1] = global_index;
		}
		if (!AmMaster)
		{
			std::vector<double> box_location(4,0);
			box_location[0] = 1.5 * (vector_comp[0]-1);
			box_location[1] = 1.5 * (vector_comp[0]);
			box_location[2] = 1.5 * y_comp;
			box_location[3] = 1.5 * (y_comp + 1);

			Box new_box(box_location);
			mHaloBoxesLeft.push_back(new_box);

			// Update the map between global and local box indices.
			unsigned global_index = vector_comp[0] - 1 + 200 * y_comp;
			mHaloBoxLeftGlobalToLocalIndexMapping[global_index] = mHaloBoxesLeft.size()-1;
			mHaloBoxLeftLocalToGlobalIndexMapping[mHaloBoxesLeft.size()-1] = global_index;
		}
	}	
	
	// Save the indices of boxes on this process which are halos.
	std::vector<unsigned> mHalosRight;
	std::vector<unsigned> mHalosLeft;
	for (unsigned y_comp = 0; y_comp<200; y_comp++)
	{
		if (!AmMaster)
		{
			mHalosLeft.push_back(vector_comp[0]+y_comp*200);
		}
		if (!AmTopMost)
		{
			mHalosRight.push_back(vector_comp[vector_comp.size()-1]+y_comp*200);
		}	
	}

	/*
	 * Set up initial Particles
	 */
	
	unsigned TotalParticles = 100000;
	std::vector<Particle*> InitialParticles;

	for (unsigned i = 0; i < TotalParticles; i++)
	{
		// Choose a random particle locations.
		std::vector<double> new_location(2,0);
		new_location[1] = DomainSize*rand()/RAND_MAX;
		new_location[0] = DomainSize*rand()/RAND_MAX;
		
		Particle new_Particle(i, new_location);
		InitialParticles.push_back(new Particle(i, new_location));
	}
	
	// Vector to store particles on each process.
	std::vector< Particle* > Particles;
	for (unsigned i = 0; i<TotalParticles; i++)
	{
		if (IsParticleOwned(InitialParticles[i], mBoxGlobalToLocalIndexMapping))
		{
			Particles.push_back(InitialParticles[i]);
		}
	}

	/******************************************
	* INITIALIZE MEMORY FOR MPI COMMUNICATIONS
	******************************************/

	// Initialise Halos and a list of those particles that are included for reference.

	std::vector< Particle* > HaloParticlesLeft;
	std::vector< Particle* > HaloParticlesRight;
	std::vector< Particle* > HaloParticlesReceivedRight;
	std::vector< Particle* > HaloParticlesReceivedLeft;

	std::vector<Particle*> LocalParticles;
	
	std::vector<Particle*> ParticlesToMoveRight,ParticlesToMoveLeft;
	std::vector<double> Space(2,0);
	

	// Initialise manual sending and receiving buffer to coniguize std::vector type matrices
	unsigned SendMessageSizeRight = 0, ReceiveMessageSizeRight = 0, SendMessageSizeLeft = 0, ReceiveMessageSizeLeft = 0;
	std::vector<double> SendRightBuffer, ReceiveRightBuffer, SendLeftBuffer, ReceiveLeftBuffer;
	std::vector<unsigned> SendIndexRightBuffer, ReceiveIndexRightBuffer, SendIndexLeftBuffer, ReceiveIndexLeftBuffer;

	/***********************************
	* Time Loop Simulation Begins Here
	************************************/
	//Some force related variables
	double SpringConstant = 0.1;
	double SpringRestLength = 0.9;
	
	
	for (unsigned t = 0; t < TimeSteps; t++)
	{

		// Set up total force map for process particles and initialize to zero.
		std::map<Particle*, std::vector<double> > Force;
		for (unsigned i = 0; i<Particles.size(); i++)
		{
			std::vector<double> zero_force(2,0);
			Force[Particles[i]] = zero_force;
		}



		// Put Particles into boxes.
		for (unsigned i = 0; i<Particles.size(); i++)
		{
			unsigned containing_box = CalculateContainingBox(Particles[i]->rGetLocation());
			// Find out which local box this is
			unsigned local_index = mBoxGlobalToLocalIndexMapping[containing_box];
			mBoxes[local_index].AddParticle(Particles[i]);
		}

		// Populate Halos if using more than one process
		
		// Clear memory
		HaloParticlesRight.clear();
		HaloParticlesLeft.clear();

		if ( numprocs > 1 )
		{
			for (unsigned i = 0; i<mHalosRight.size(); i++)
			{
				for (std::set<Particle*>::iterator iter = mBoxes[mBoxGlobalToLocalIndexMapping[mHalosRight[i]]].rGetParticlesContained().begin();
					iter !=  mBoxes[mBoxGlobalToLocalIndexMapping[mHalosRight[i]]].rGetParticlesContained().end();
					++iter)
				{
					HaloParticlesRight.push_back(*iter);
				}
			}
			for (unsigned i = 0; i<mHalosLeft.size(); i++)
			{
				for (std::set<Particle*>::iterator iter = mBoxes[mBoxGlobalToLocalIndexMapping[mHalosLeft[i]]].rGetParticlesContained().begin();
					iter !=  mBoxes[mBoxGlobalToLocalIndexMapping[mHalosLeft[i]]].rGetParticlesContained().end();
					++iter)
				{
					HaloParticlesLeft.push_back(*iter);
				}
			} 	
		}	

		SendMessageSizeRight = 2*HaloParticlesRight.size();	
		SendMessageSizeLeft = 2*HaloParticlesLeft.size();

		begin = MPI_Wtime();

		MPI_Sendrecv(	&SendMessageSizeRight,
				1, 
				MPI_UNSIGNED, 
				right, 
				t, 
				&ReceiveMessageSizeLeft, 
				1, 
				MPI_UNSIGNED, 
				left, 
				t, 
				MPI_COMM_WORLD, 
				&status);

		MPI_Sendrecv(	&SendMessageSizeLeft, 
				1, 
				MPI_UNSIGNED, 
				left, 
				t, 
				&ReceiveMessageSizeRight, 
				1, 
				MPI_UNSIGNED, 
				right, 
				t, 
				MPI_COMM_WORLD, 
				&status);

		end = MPI_Wtime();
		mpi_time+= (end-begin);


		// Clear recv buffer and resize.
		ReceiveLeftBuffer.clear();
		ReceiveRightBuffer.clear();
		ReceiveLeftBuffer.resize(ReceiveMessageSizeLeft);
		ReceiveRightBuffer.resize(ReceiveMessageSizeRight);

		ReceiveIndexLeftBuffer.clear();
		ReceiveIndexRightBuffer.clear();
		ReceiveIndexLeftBuffer.resize(ReceiveMessageSizeLeft/2);
		ReceiveIndexRightBuffer.resize(ReceiveMessageSizeRight/2);

		// Buffer halo data into contiguous form
		SendLeftBuffer.clear();
		SendRightBuffer.clear();
		SendRightBuffer.resize(SendMessageSizeRight);
		SendLeftBuffer.resize(SendMessageSizeLeft);

		SendIndexLeftBuffer.clear();
		SendIndexRightBuffer.clear();
		SendIndexRightBuffer.resize(SendMessageSizeRight/2);
		SendIndexLeftBuffer.resize(SendMessageSizeLeft/2);

		// Store Halo data as for example (x1,y1,x2,y2,...) for 2D 
		for (unsigned i = 0; i<HaloParticlesRight.size(); i++)
		{
			SendIndexRightBuffer[i] = HaloParticlesRight[i]->GetIndex();
			for (unsigned d = 0;d<2;d++)
			{
				SendRightBuffer[2*i+d] = HaloParticlesRight[i]->rGetLocation()[d];
			}

		}
		for (unsigned i = 0; i<HaloParticlesLeft.size(); i++)
		{
			SendIndexLeftBuffer[i] = HaloParticlesLeft[i]->GetIndex();
			for (unsigned d = 0;d<2;d++)
			{
				SendLeftBuffer[2*i+d] = HaloParticlesLeft[i]->rGetLocation()[d];
			}
		}


		begin = MPI_Wtime();

		// Send and receive halo data
		MPI_Sendrecv(	&SendRightBuffer[0], 
				SendMessageSizeRight, 
				MPI_DOUBLE, 
				right, 
				t, 
				&ReceiveLeftBuffer[0],
				ReceiveMessageSizeLeft, 
				MPI_DOUBLE, 
				left, 
				t, 
				MPI_COMM_WORLD, 
				&status);

		MPI_Sendrecv(	&SendIndexRightBuffer[0], 
				SendMessageSizeRight/2, 
				MPI_UNSIGNED, 
				right, 
				t, 
				&ReceiveIndexLeftBuffer[0],
				ReceiveMessageSizeLeft/2, 
				MPI_UNSIGNED, 
				left,
				t, 
				MPI_COMM_WORLD, 
				&status);
		
		MPI_Sendrecv(	&SendLeftBuffer[0], 
				SendMessageSizeLeft, 
				MPI_DOUBLE, 
				left, 
				t, 
				&ReceiveRightBuffer[0],
				ReceiveMessageSizeRight,
				MPI_DOUBLE, 
				right,
				t, 
				MPI_COMM_WORLD, 
				&status);
	
		MPI_Sendrecv(	&SendIndexLeftBuffer[0], 	
				SendMessageSizeLeft/2, 
				MPI_UNSIGNED, 
				left, 
				t, 
				&ReceiveIndexRightBuffer[0],
				ReceiveMessageSizeRight/2, 
				MPI_UNSIGNED, 
				right,
				t, 
				MPI_COMM_WORLD, 
				&status);

		end = MPI_Wtime();
		mpi_time+= (end-begin);

		// Make halo particles
		HaloParticlesReceivedLeft.clear();
		HaloParticlesReceivedRight.clear();

		// Clear halo boxes
		for (unsigned i = 0; i<mHaloBoxesLeft.size(); i++)
		{
			mHaloBoxesLeft[i].Clear();
		}
		for (unsigned i = 0; i<mHaloBoxesRight.size(); i++)
		{
			mHaloBoxesRight[i].Clear();
		}

		for (unsigned i = 0; i<ReceiveMessageSizeLeft/2; i++)
		{
			std::vector<double> location(2,0);
			for (unsigned d = 0;d<2;d++)
			{
				location[d] = ReceiveLeftBuffer[2*i+d];
			}
			Particle* p_Particle = new Particle(ReceiveIndexLeftBuffer[i] , location);
			HaloParticlesReceivedLeft.push_back(p_Particle);

			// Add it to the local halo box.
			unsigned local_index = mHaloBoxLeftGlobalToLocalIndexMapping[CalculateContainingBox(location)];
			mHaloBoxesLeft[local_index].AddParticle(p_Particle);
		}
		for (unsigned i = 0; i<ReceiveMessageSizeRight/2; i++)
		{
			std::vector<double> location(2,0);
			for (unsigned d = 0; d<2; d++)
			{
				location[d] = ReceiveRightBuffer[2*i+d];
			}
			Particle* p_Particle = new Particle(ReceiveIndexRightBuffer[i] , location);
			HaloParticlesReceivedRight.push_back(p_Particle);

			// Add it to the local halo box.
			unsigned local_index = mHaloBoxRightGlobalToLocalIndexMapping[CalculateContainingBox(location)];
			mHaloBoxesRight[local_index].AddParticle(p_Particle);
		}

		// Calculate repulsion force on all local particles.
		std::map<Particle*, std::vector<double> > repulsion_force;
		double START = MPI_Wtime();
		// For each of the particles find neighbouring particles and add a contribution to the force on that particle
		for (unsigned i = 0; i<Particles.size(); i++)
		{
			std::vector<double> location = Particles[i]->rGetLocation(); 

			// Get containing box.
			unsigned containing_box = CalculateContainingBox(Particles[i]->rGetLocation());
			
			// Get local index of that box
			unsigned local_index = mBoxGlobalToLocalIndexMapping[containing_box];

			// Loop over the local boxes
			for (std::set<unsigned>::iterator local_box_iter = mLocalBoxes[local_index].begin();
				local_box_iter !=  mLocalBoxes[local_index].end();
				++local_box_iter)
			{
				// Find where this box lives (local or halo left / right).
				bool is_left_halo = (mHaloBoxLeftGlobalToLocalIndexMapping.find(*local_box_iter) !=  mHaloBoxLeftGlobalToLocalIndexMapping.end() );
 				bool is_right_halo = (mHaloBoxRightGlobalToLocalIndexMapping.find(*local_box_iter) !=  mHaloBoxRightGlobalToLocalIndexMapping.end() );

				// Get the set of Particles contained
				std::set<Particle*> Particles_contained;
				if (is_left_halo)
				{
					Particles_contained = mHaloBoxesLeft[mHaloBoxLeftGlobalToLocalIndexMapping[*local_box_iter]].rGetParticlesContained();
				}
				else if (is_right_halo)
				{
					Particles_contained = mHaloBoxesRight[mHaloBoxRightGlobalToLocalIndexMapping[*local_box_iter]].rGetParticlesContained();
				}
				else
				{
					Particles_contained = mBoxes[mBoxGlobalToLocalIndexMapping[*local_box_iter]].rGetParticlesContained();
				}

				// Loop over each of these Particles
				for (std::set<Particle*>::iterator Particle_iter = Particles_contained.begin();
					Particle_iter !=  Particles_contained.end();
					++Particle_iter)
				{
					// Add force contribution
					std::vector<double> location1 = (*Particle_iter)->rGetLocation();
					std::vector<double> separation(2,0);
					separation[0] = location1[0] - location[0];
					separation[1] = location1[1] - location[1];

					double distance = sqrt( pow((location[0] - location1[0]),2) + pow((location[1] - location1[1]),2));
					
					// Normalize the separation
					separation[0] = separation[0] / distance;
					separation[1] = separation[1] / distance;
					
					if (distance > 0 && distance < 1.5)
					{
						Force[Particles[i]][0] +=  (-1)*SpringConstant*(SpringRestLength-distance)*separation[0];
						Force[Particles[i]][1] +=  (-1)*SpringConstant*(SpringRestLength-distance)*separation[1]; 
					}
				}

			}
		}

		double END = MPI_Wtime();
		time_calculating_pairs+= (END-START);

		
		// Calculate random motion force
		std::vector< std::vector<double> > random_force(Particles.size(), std::vector<double>(2,0));
		// Set diffusion constant
		double D = 0.006;
		random_force = CalculateRandomForce(Particles, DampingConstant, DeltaT, D);

		// Add random force
		for (unsigned i = 0; i<Particles.size(); i++)
		{
			for (unsigned d = 0; d < 2; d++)
			{
				Force[Particles[i]][d]+= random_force[i][d];
			}
		}

		// Move particles (see DOI: 10.1046/j.0960-7722.2001.00216.x for details)

		double constant = DeltaT/DampingConstant;
	
		// Clear the containers of particles that need to move.
		ParticlesToMoveLeft.clear();
		ParticlesToMoveRight.clear();
		for (unsigned i = 0; i < Particles.size(); i++)
		{
			double new_x_location = Particles[i]->rGetLocation()[0] + constant * Force[Particles[i]][0];
			double new_y_location = Particles[i]->rGetLocation()[1] + constant * Force[Particles[i]][1];

			// Apply reflective boundary conditions.
			if (new_x_location < 0.0)			
			{
				new_x_location = 2 * 0.0 - new_x_location;
				if (new_y_location < 0.0)
				{
					new_y_location = 2 * 0.0 - new_y_location; 
				}
				else if (new_y_location > DomainSize)
				{
					new_y_location = 2 * DomainSize - new_y_location; 
				}
				else
				{
					// Do nothing to y-vector.
				}	
			}
			else if (new_x_location > DomainSize)
			{
				new_x_location = 2 * DomainSize - new_x_location;
				if (new_y_location < 0.0)
				{
					new_y_location = 2 * 0.0 - new_y_location; 
				}
				else if (new_y_location > DomainSize)
				{
					new_y_location = 2 * DomainSize - new_y_location; 
				}
				else
				{
					// Do nothing to y-vector.
				}
			}
			else
			{
				if (new_y_location < 0.0)
				{
					new_y_location = 2 * 0.0 - new_y_location; 
				}
				else if (new_y_location > DomainSize)
				{
					new_y_location = 2 * DomainSize - new_y_location; 
				}
				else
				{
					// Do nothing to y-vector.
				}			
			}

			Particles[i]->rGetLocation()[0] = new_x_location;
			Particles[i]->rGetLocation()[1] = new_y_location;

			// Check particles still lie in area for the given process, and flag any 
			// Particles that have moved off the process to be sent.

			// Only move particles across processes if np>1
			if (numprocs > 1)
			{
				if (Particles[i]->rGetLocation()[0]<RegionBound[0])
				{
					ParticlesToMoveLeft.push_back(Particles[i]);
				}
				if (Particles[i]->rGetLocation()[0]>RegionBound[1] || Particles[i]->rGetLocation()[0] == RegionBound[1])
				{
					ParticlesToMoveRight.push_back(Particles[i]);
				}
			}
		}
		

		// Send particle position data to process right and left
		SendMessageSizeLeft  = 2 * ParticlesToMoveLeft.size();
		SendMessageSizeRight = 2 * ParticlesToMoveRight.size();
		
		SendRightBuffer.clear();
		SendLeftBuffer.clear();
		SendRightBuffer.resize(SendMessageSizeRight);
		SendLeftBuffer.resize(SendMessageSizeLeft);

		SendIndexRightBuffer.clear();
		SendIndexLeftBuffer.clear();
		SendIndexRightBuffer.resize(SendMessageSizeRight/2);
		SendIndexLeftBuffer.resize(SendMessageSizeLeft/2);

		for (unsigned i = 0; i<ParticlesToMoveLeft.size(); i++)
		{
			SendIndexLeftBuffer[i] = ParticlesToMoveLeft[i]->GetIndex();
			for (unsigned d = 0; d < 2; d++)
			{
				SendLeftBuffer[2*i+d] = ParticlesToMoveLeft[i]->rGetLocation()[d];
			}
		}

		for (unsigned i = 0; i<ParticlesToMoveRight.size(); i++)
		{
			SendIndexRightBuffer[i] = ParticlesToMoveRight[i]->GetIndex();
			for (unsigned d = 0; d < 2; d++)
			{
				SendRightBuffer[2*i+d] = ParticlesToMoveRight[i]->rGetLocation()[d];
			}
		}

		begin = MPI_Wtime();

		MPI_Sendrecv(&SendMessageSizeLeft, 1, MPI_UNSIGNED, left, t, &ReceiveMessageSizeRight, 1, MPI_UNSIGNED, right, t, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&SendMessageSizeRight, 1, MPI_UNSIGNED, right, t, &ReceiveMessageSizeLeft, 1, MPI_UNSIGNED, left, t, MPI_COMM_WORLD, &status);

		end = MPI_Wtime();
		mpi_time+= (end-begin);

		ReceiveLeftBuffer.clear();
		ReceiveRightBuffer.clear();
		ReceiveLeftBuffer.resize(ReceiveMessageSizeLeft);
		ReceiveLeftBuffer.resize(ReceiveMessageSizeRight);

		ReceiveIndexLeftBuffer.clear();
		ReceiveIndexRightBuffer.clear();
		ReceiveIndexLeftBuffer.resize(ReceiveMessageSizeLeft/2);
		ReceiveIndexRightBuffer.resize(ReceiveMessageSizeRight/2);

		begin = MPI_Wtime();

		MPI_Sendrecv(&SendRightBuffer[0], SendMessageSizeRight , MPI_DOUBLE, right, t, &ReceiveLeftBuffer[0], ReceiveMessageSizeLeft , MPI_DOUBLE, left, t, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&SendIndexRightBuffer[0], SendMessageSizeRight/2 , MPI_UNSIGNED, right, t, &ReceiveIndexLeftBuffer[0], ReceiveMessageSizeLeft/2 , MPI_UNSIGNED, left, t, MPI_COMM_WORLD, &status);
		
		MPI_Sendrecv(&SendLeftBuffer[0], SendMessageSizeLeft , MPI_DOUBLE, left, t, &ReceiveRightBuffer[0] ,ReceiveMessageSizeRight , MPI_DOUBLE, right, t, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&SendIndexLeftBuffer[0], SendMessageSizeLeft/2 , MPI_UNSIGNED, left, t, &ReceiveIndexRightBuffer[0] ,ReceiveMessageSizeRight/2 , MPI_UNSIGNED, right, t, MPI_COMM_WORLD, &status);

		end = MPI_Wtime();
		mpi_time+= (end-begin);

		// Delete moved particles from current process.
		
		// Make a copy of particles
		std::vector<Particle*> particle_copy(Particles);
		
		// Clear Particles
		Particles.clear();
		assert(Particles.size() ==  0);
		Particles.resize(particle_copy.size() - ParticlesToMoveRight.size() - ParticlesToMoveLeft.size() + ReceiveMessageSizeRight/2 + ReceiveMessageSizeLeft/2);

		// Replace un-deleted particles.
		unsigned spacer = 0;
		for (unsigned i = 0; i<particle_copy.size(); i++)
		{
			bool is_moved_left = false;
			bool is_moved_right = false;

			// Efficiency of this assumes that very few particles move across the boundary at a time-step
			for ( unsigned j = 0; j<ParticlesToMoveRight.size(); j++)
			{
				if ((particle_copy[i] ==  ParticlesToMoveRight[j]))
				{
					is_moved_right = true;
				}
			}
			for ( unsigned j = 0; j<ParticlesToMoveLeft.size(); j++)
			{
				if ((particle_copy[i] ==  ParticlesToMoveLeft[j]))
				{
					is_moved_left = true;
				}
			}

			if (!is_moved_left && !is_moved_right)
			{
				Particles[spacer] = particle_copy[i];
				spacer++;
			}
		}

		// Create new Particles from received data
		for (unsigned i = 0; i<ReceiveMessageSizeLeft/2; i++)
		{ 
			unsigned index = ReceiveIndexLeftBuffer[i];
			std::vector<double> location(2,0);
			location[0] = ReceiveLeftBuffer[2*i];
			location[1] = ReceiveLeftBuffer[2*i+1];
			Particle* p_new_Particle = new Particle(index, location);
			Particles[spacer] = p_new_Particle;
			spacer++;
		}
		for (unsigned i = 0; i<ReceiveMessageSizeRight/2; i++)
		{
			unsigned index = ReceiveIndexRightBuffer[i];
			std::vector<double> location(2,0);
			location[0] = ReceiveRightBuffer[2*i];
			location[1] = ReceiveRightBuffer[2*i+1];
			Particle* p_new_Particle = new Particle(index, location);
			Particles[spacer] = p_new_Particle;
			spacer++;
		}

		// Make sure we have filled up particles
		assert(spacer ==  Particles.size());

		for (unsigned i = 0; i<Particles.size(); i++)
		{
			assert(Particles[i]->rGetLocation()[0] > -1.0 &&  Particles[i]->rGetLocation()[0]< 301.0);
			assert(Particles[i]->rGetLocation()[1] > -1.0 &&  Particles[i]->rGetLocation()[1]< 301.0);
		}

		// Check here that particle numbers are conserved.
		unsigned num_local_particles = Particles.size();
		unsigned total_particles;
		MPI_Allreduce(&num_local_particles ,&total_particles, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
		
		assert(total_particles = TotalParticles);
		
		// If output is swtiched on (see line 387) then collect data using MPI reduction.
		if (WriteOutputToFiles)
		{		
			// Buffer all data on each process for outputing.
			double ParticlesBuffer[4*Particles.size()];

			for (unsigned i = 0; i<Particles.size(); i++)
			{
				for (unsigned d = 0;d<2;d++)
				{
					ParticlesBuffer[4*i+d] = Particles[i]->rGetLocation()[d];
				}
			
				// Tag with process id
				ParticlesBuffer[4*i+2] = Particles[i]->GetIndex();
				ParticlesBuffer[4*i+3] = pid;			
			}


			unsigned num_particles = Particles.size();

			// Send data to process 0.
			begin = MPI_Wtime();

			MPI_Send(&num_particles,1,MPI_UNSIGNED,0,t,MPI_COMM_WORLD);
			MPI_Send(ParticlesBuffer, 4*num_particles, MPI_DOUBLE, 0,t, MPI_COMM_WORLD);
		
			end = MPI_Wtime();
			mpi_time+= (end-begin);

			// Process 0 receieves data from each process and stores it
			if (pid == 0)
			{		
				// Memory to store all particles

				std::vector< std::vector<double> > AllParticles(TotalParticles, std::vector<double>(4,0));

				double ParticleDataBuffer[4*TotalParticles];
		
				unsigned space_used = 0;
				for (unsigned i = 0;i<numprocs;i++)
				{	
					begin = MPI_Wtime();
					unsigned temp_num_particles;
					MPI_Recv(&temp_num_particles,1,MPI_UNSIGNED,i,t,MPI_COMM_WORLD,&status);
					MPI_Recv(&ParticleDataBuffer[space_used],4*temp_num_particles,MPI_DOUBLE,i,t,MPI_COMM_WORLD,&status);

					end = MPI_Wtime();
					mpi_time+= (end-begin);
	
					space_used+= 4*temp_num_particles;
				}

				// Unpack data
				for (unsigned i = 0;i<AllParticles.size();i++)
				{
					AllParticles[i][0] = ParticleDataBuffer[4*i];
					AllParticles[i][1] = ParticleDataBuffer[4*i+1];
					AllParticles[i][2] = ParticleDataBuffer[4*i+2];
					AllParticles[i][3] = ParticleDataBuffer[4*i+3];
				}

				// pid 0 to output data
		
				file << t << " ";
				for (unsigned i = 0;i<AllParticles.size();i++)
				{
					file << " " << AllParticles[i][0] << " " << AllParticles[i][1] << " " << AllParticles[i][2] << " " << AllParticles[i][3];
			
				}
				file << "\n";
		
			}
		

			std::cout << t << "...";
		}
	}

/*	Load balance information commented for easier post-processing of scaling results.
//Output total times and times on each process.
std::cout << "Time on process " << pid << " was " << time << "\n";
std::cout << "Time in pairs on process " << pid << " was " << time_calculating_pairs << "\n";
std::cout << "Time on MPI on process " << pid << " was " << mpi_time << "\n";
*/

// MPI Timing
double finish = MPI_Wtime();
double time = finish-start;
double totaltime, totalmpi_time, pairs_total_time;


MPI_Allreduce(&time,&totaltime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
MPI_Allreduce(&mpi_time,&totalmpi_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(&time_calculating_pairs,&pairs_total_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

if (pid == 0)
{
	std::cout << "Total time was " << totaltime <<"\n";
	std::cout << "Total MPI time was " << totalmpi_time << "\n";
}	


MPI::Finalize();
return 0;
}


