// Symmetric Toth XMPF written with MMSP classes and no optimizations

#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE

#ifdef BGQ
#include<mpi.h>
#endif
#include<iomanip>
#include<vector>
#include<cmath>
#include<ctime>
#include<cstdio>
#include<cstdlib>
#include<cassert>
#include"MMSP.hpp"

typedef float phi_type; 		// Note : This is used extensively in sp-initialize.cpp
typedef MMSP::sparse<phi_type> store_type; // Note : This is used extensively in sp-initialize.cpp

//included after typedefs for proper variable types in the functions of the included files
#include"sp-graingrowth.hpp" 
#include"sp-initialize.cpp"

//Function declarations 
std::vector<std::string> split(const std::string &text, char sep);
void print_progress(const int step, const int steps, const int iterations);
template <int dim> void makeSplit(MMSP::grid<dim, store_type >& grid);

namespace MMSP {

void generate(int dim, char* filename) {
//	Use case for this code is as follows: generate() will create MMSP grid files which can be used to store field data.
//	The type of grid data to be loaded is selected by the name of the filename that you declare, for example: 
//		 planar.000000.dat to generate a two-grain, planar interface
//		 circle.000000.dat to generate a radially-symmetric grain-in-grain structure
//		 filexXxYxZxN.000000.dat to read an appropriatelly generated sharp-interface input file 
//			(dimensions X, Y, Z and with N grains eg filex100x100x15.000.dat for a grid of 100x100 with 15 grains)

	std::string search_name(filename);
	// Check if the filename has each of the following phrases in them. string search function return npos if no matches are found. npos here refers to the end of file. 
	bool planar = (search_name.find("planar") != std::string::npos);
	bool circle = (search_name.find("circle") != std::string::npos);
	bool generated = (search_name.find("file") != std::string::npos);
	
	int nx = 100;    //Default values of nx and ny for the two dimensions.
	int ny = 100;
	int num_grains = 2;
	std::vector<int> dimensions;	//Each component of dimensions stores the number of grid points in that dimension.
	if (generated){// determining the dimensions based off filename data
		std::vector<std::string> splits = split(search_name, '.');
		std::string name_root = splits[0]; //"filex#####x#####
		std::vector<std::string> metadata = split(name_root, 'x');
		for (int i = 0; i < dim; i++){
			dimensions.push_back(std::atoi(metadata[i+1].c_str()));     //atoi interprets a string of numbers as integer values. c_str converts the std::string of C++ into a C string
		}
		num_grains = std::atoi(metadata[dim+1].c_str());
		nx = dimensions[0];
		ny = dimensions[1];
	}
	else{
		for (int i = 0; i < dim; i++){
			dimensions.push_back(0);   			
		}
	}

	if(planar){
		MMSP::grid<2,store_type > grid (2,0,nx,0,ny);    //Define the MMSP grid. The first input parameter refers to the number of grains.
	
		MMSP::dx(grid, 0) = 1.0; //Lx/nx ; Thus, Lx = Ly = 100.0 ;
		MMSP::dx(grid, 1) = 1.0; //Ly/ny;
		std::cout << "Creating 2-grain planar interface."<<std::endl;
		
		makeSplit<2>(grid);
		
		std::cout << "Saving..." << std::endl; 
		output(grid, filename); //write out initialized grid 
		std::cout << "Grid saved as " << filename << std::endl;
	} else if(circle){
		MMSP::grid<2,store_type > grid (2,0,nx,0,ny);  //Define the MMSP grid. The first input parameter refers to the number of grains.
	
		MMSP::dx(grid, 0) = 1.0;
		MMSP::dx(grid, 1) = 1.0;
		std::cout << "Creating circular grain-in-grain test grid."<<std::endl;
		
		float radius = 20;   //Define the radius by default.
		makeCenterCircle<2>(grid, radius);
		
		std::cout << "Saving..." << std::endl;
		output(grid, filename); //write out initialized grid
		std::cout << "Grid saved as " << filename << std::endl;
	} else if(generated){
		MMSP::grid<2,store_type > grid (num_grains,0,nx,0,ny);  //Define the MMSP grid. The first input parameter refers to the number of grains.
	
		MMSP::dx(grid, 0) = 1.0;//Lx/nx;
		MMSP::dx(grid, 1) = 1.0;//Ly/ny;

		std::cout << "Reading in generated microstructure."<<std::endl;
		
		treadin<2,store_type>(grid);
		
		std::cout << "Saving..." << std::endl;
		output(grid, filename); //write out initialized grid
		std::cout << "Grid saved as " << filename << std::endl;
	}
	else {
		std::cerr<<"Error: No initialization condition selected."<<std::endl;
		exit(1);
	}
	
}

template <int dim> void update(MMSP::grid<dim, store_type >& grid, int steps) {    //computation of the model 
	float dx = MMSP::dx(grid, 0);   // functor
	int id=0;
	#ifdef MPI_VERSION
 	id=MPI::COMM_WORLD.Get_rank();     
	#endif
	
	//thresh is the minimum phi value imposed to minimize mathematical noise in fields
	phi_type thresh = 1.0e-3;
	
	phi_type width = 1.0;				//width of stable interface
	phi_type M = 1.0e-3;				//mobility of field motion
	phi_type dt = dx / (20.0 * M);		//timestep calculated as function of resolution, 10x less than Courant maximum
	phi_type gamma = 1.0/3.0;			//energy contribution from gradient component
	phi_type w = 3.0 * (gamma/width);	//energy contribution from double-well component
	phi_type eps_sq = 3.0 * width * gamma;
	
// Spatially-varying simulation parameters calculated
// Initialize vairiables
	const int num_grains =  fields(grid); 	//functor
	
	static int iterations = 1;
	
	for (int step = 0; step < steps; step++) {
		if(id == 0)	print_progress(step, steps, iterations);
		ghostswap(grid);    //exchange of information 
		
		MMSP::grid<dim, store_type > update(grid);	// ?? Potential for improvement : lots of duplications which deter easy understanding of MMSP

		for (int n = 0; n < nodes(grid); n++){		
			int i = n;
			vector<int> x = position(grid, n);
			
			// determine nonzero fields within
			// the neighborhood of this node
			// (2 adjacent voxels along each cardinal direction)
			sparse<int> s;
			for (int j = 0; j < dim; j++){
				for (int k = -1; k <= 1; k++) {
				  x[j] += k;
				  for (int h = 0; h < length(grid(x)); h++) {
				    int index = MMSP::index(grid(x), h);
				    set(s, index) = 1;
				  }
				  x[j] -= k;
				}
			}
			phi_type S = phi_type(length(s));

			// if only one field is nonzero,
			// then copy this node to update
			if (S < 2.0){
				update(n) = grid(n);
			} else {
				// compute laplacian of each field
				sparse<phi_type> lap = laplacian(grid, n);
				
				sparse<phi_type> dFdp;

				for (int h = 0; h < length(s); h++) {
					// compute variational derivatives
					int hindex = MMSP::index(s, h);
					phi_type N = num_grains;
					
					phi_type phi_sq = 0.0; // phi_sq is the sum of the squares of the phi field
					phi_type dFall = 0.0;
					
					// Compute phi_sq value by taking the sum of all (defined) phi values
					for (int j = 0; j < length(s); j++) {
						int jindex = MMSP::index(s, j);
						phi_sq += grid(n)[jindex]*grid(n)[jindex];
					}
					
					// Compute the dFall value by calculating all (defined) dFdp values and summing them
					for (int j = 0; j < length(s); j++) {
						int jindex = MMSP::index(s, j);
						set(dFdp, jindex) = w * grid(n)[jindex] * (phi_sq - grid(n)[jindex]) - eps_sq * (lap[jindex]);
						dFall += dFdp[jindex];
					}
					
					store_type dpdt;
					set(dpdt, hindex) = - (M/N) * (N*dFdp[hindex] - dFall);
					phi_type value = grid(n)[hindex] + dt * dpdt[hindex];
					if (value > 1.0) value = 1.0;
					if (value < -0.001) value = 0.0;
					if (value > thresh) set(update(n), hindex) = value;
					else if (grid(n)[hindex] != 0.0) set(update(n), hindex) = 0.0;
				}//calculate interactions between interacting fields
			}//perform calculations on non-zero fields
		} //loop over nodes		
		swap(grid, update);
	} //loop over steps
	ghostswap(grid);
	++iterations;
}

template <class T> std::ostream& operator<<(std::ostream& o, sparse<T>& s) {
	o<<"    Index  Value\n";
	for (int i=0; i<length(s); ++i) {
		int index = MMSP::index(s, i);
		o<<"    "<<std::setw(5)<<std::right<<index<<"  "<<s[index]<<'\n';
	}
	return o;
}

} // namespace MMSP

void print_progress(const int step, const int steps, const int iterations) {
	char* timestring;
	static unsigned long tstart;
	struct tm* timeinfo;

	if (step==0) {
		tstart = time(NULL);
		std::time_t rawtime;
		std::time( &rawtime );
		timeinfo = std::localtime( &rawtime );
		timestring = std::asctime(timeinfo);
		timestring[std::strlen(timestring)-1] = '\0';
		std::cout<<"Pass "<<std::setw(3)<<std::right<<iterations<<": "<<timestring<<" ["<<std::flush;
	} else if (step==steps-1) {
		unsigned long deltat = time(NULL)-tstart;
		std::cout<<"•] "<<std::setw(2)<<std::right<<deltat/3600<<"h:"
										<<std::setw(2)<<std::right<<(deltat%3600)/60<<"m:"
										<<std::setw(2)<<std::right<<deltat%60<<"s"
										<<" (File "<<std::setw(5)<<std::right<<iterations*steps<<")."<<std::endl;
	} else if ((20 * step) % steps == 0) std::cout<<"• "<<std::flush;
}

std::vector<std::string> split(const std::string &text, char sep) {
  std::vector<std::string> tokens;
  std::size_t start = 0, end = 0;
  while ((end = text.find(sep, start)) != std::string::npos) {
    tokens.push_back(text.substr(start, end - start));
    start = end + 1;
  }
  tokens.push_back(text.substr(start));
  return tokens;
}

#endif

#include"MMSP.main.hpp"
