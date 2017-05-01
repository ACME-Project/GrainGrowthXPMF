// Tool take 2

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>

#include <zlib.h>
#include <IL/il.h>
#include <IL/ilu.h>
#include <IL/ilut.h>
#include "devil_cpp_wrapper.hpp"

#include "MMSP.hpp"

void writePNG(const int w, const int h, const int bpp, unsigned char* imData, const char* filename, bool outline);
std::vector<std::string> split(const std::string &text, char sep);

int main(int argc, char* argv[]) {

	std::ifstream input(argv[1]);
	if (!input) {
		std::cerr<<"File input error: could not open "<<argv[1]<<".\n\n";
		exit(-1);
	}
	
  // determining the step number for the current file
	const std::string name = std::string(argv[1]);
	std::vector<std::string> splits = split(name, '.');
	std::string name_root = splits[0];
	std::string number = splits[1];
	std::string output_dir = std::string(argv[2]);
	int step_num = std::atoi(number.c_str());
	
// file header metadata scraping "Describe"
	std::string type;
	getline(input,type,'\n');
	bool sparse_type = (type.find("sparse") != std::string::npos);
	int dim;
	input>>dim;
	int fields;
	input >> fields;
	int x0[3] = {0, 0, 0};
	int x1[3] = {0, 0, 0};
	for (int i = 0; i < dim; i++)
	input >> x0[i] >> x1[i];
	float dx[3] = {1.0, 1.0, 1.0};
	for (int i = 0; i < dim; i++)
	input >> dx[i];
	int width=x1[0]-x0[0];
	int height=x1[1]-x0[1];
	
	
	unsigned char* buffer = new unsigned char[(height*width)];
	std::vector<int> grain_IDs; // list of grain id #s
	std::vector<float> grain_areas; //list of grain areas
	float dA = 1.0;
	float threshold = 1.0;

	if(sparse_type){
		MMSP::grid<2,MMSP::sparse<float> > grid(argv[1]);
// file data collection "Detect"
		for (int n=0; n<nodes(grid); n++) {
			MMSP::vector<int> x = position(grid, n);
		
			MMSP::sparse<float> this_node = grid(n);
			float sum = 0.0;
			for (int l = 0; l < MMSP::length(this_node); l++){
				int dex = MMSP::index(this_node, l);
				float phi = this_node[dex];
				sum += phi*phi;
		
				bool found = false;
				for(int id_pos = 0; id_pos < grain_IDs.size(); id_pos++){
					if (dex == grain_IDs[id_pos]){
						grain_areas[id_pos] += phi*dA;
						found = true;
					}
				}
				if (!found){
					grain_IDs.push_back(dex);
					grain_areas.push_back(phi);
				}
			
			}
			buffer[n] = 255*sqrt(sum);
	//		std::cout<<"Total: "<< sum<<std::endl;
		}	
	} else {
		MMSP::grid<2,MMSP::vector<float> > grid(argv[1]);
// file data collection "Detect"
		for (int n=0; n<nodes(grid); n++) {
			MMSP::vector<int> x = position(grid, n);
		
			MMSP::vector<float> this_node = grid(n);
			float sum = 0.0;
			for (int l = 0; l < MMSP::length(this_node); l++){
				int dex = l;
				float phi = this_node[dex];
				sum += phi*phi;
				buffer[n] = 255*sqrt(sum);
				bool found = false;
				for(int id_pos = 0; id_pos < grain_IDs.size(); id_pos++){
					if (dex == grain_IDs[id_pos]){
						grain_areas[id_pos] += phi*dA;
						found = true;
					}
				}
				if (!found){
					grain_IDs.push_back(dex);
					grain_areas.push_back(phi);
				}
			
			}
	//		std::cout<<"Total: "<< sum<<std::endl;
		}	
	}
	
	std::ostringstream outs;
	outs << output_dir<< name_root << "." << number << "_areas.csv";
	std::string oname = outs.str();
	std::ofstream outfile;
	outfile.open(oname.c_str());
	for (int grain_number = 0; grain_number<grain_IDs.size(); grain_number++){
		if (grain_areas[grain_number] > threshold) outfile <<grain_areas[grain_number] << std::endl;
	}
	outfile.close();
	std::cout<<"Grain data written to "<<oname<<"."<<std::endl; 

	std::ostringstream bs;
	bs << output_dir<< name_root << "." << number << "_grains.png";
	std::string png_name = bs.str();
	writePNG(width, height, 1, buffer, png_name.c_str(), true);
	std::cout<<"Microstructure visualized in "<<png_name<<"."<<std::endl; 
	
	return 0;
}


void writePNG(const int w, const int h, const int bpp, unsigned char* imData, const char* filename, bool outline=true)
{
  // Initialize image
 	ilInit();
 	ILenum Error;
  ILuint imageID = ilGenImage() ;
 	ilBindImage(imageID);
 	if (outline) ilTexImage(h, w, 1, bpp, IL_LUMINANCE, IL_UNSIGNED_BYTE, imData);
  	else ilTexImage(h, w, 1, bpp, IL_RGB, IL_UNSIGNED_BYTE, imData);
//  ilTexImage(w, h, 1, bpp, IL_LUMINANCE, IL_UNSIGNED_BYTE, imData); //original line from Trevor
 	Error = ilGetError();
  if (Error!=IL_NO_ERROR) std::cout<<"Error making image: "<<iluErrorString(Error)<<std::endl;
 	ilEnable(IL_FILE_OVERWRITE);
  ilSave( IL_PNG, filename) ;
 	Error = ilGetError();
  if (Error!=IL_NO_ERROR) std::cout<<"Error saving image: "<<iluErrorString(Error)<<std::endl;
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

