#include <fstream>
#include <iostream>
using namespace std;

// This file illustrates the basic reading of a Gadget snapshot (SnapFormat 2)
// We read the number of gas, DM and star particles and put their positions,
// smoothing lengths, densities and temperatures in arrays

// This method reads positions. We first read the size of the position block,
// then we need to skip 8 bytes and then we can read the positions themselves
// Finally, we skip the last 4 bytes.
void read_positions(float* positions, ifstream& ifile){
  unsigned int blocksize;
  ifile.read((char*)&blocksize, 4);
  // skip end of name block and beginning of data block
  ifile.seekg(8, ios_base::cur);
  // the number of particles is the size of the block minus the 4 bytes at the
  // start and end, divided by 3 (3 coordinates) and divided by 4 (floats)
  unsigned int numpart = (blocksize-8)/12;
  for(unsigned int i = 0; i < numpart; i++){
    // we read the 3 coordinates in 1 go
    ifile.read((char*)&positions[3*i], 12);
  }
  // come to think of it: actually we can read the entire array in 1 go...
  ifile.seekg(4, ios_base::cur);
}

// This method reads an array. The idea is exactly the same as above, apart from
// the fact that we read 1 coordinate per particle.
void read_array(float* array, ifstream& ifile){
  unsigned int blocksize;
  ifile.read((char*)&blocksize, 4);
  // skip end of name block and beginning of data block
  ifile.seekg(8, ios_base::cur);
  unsigned int numpart = (blocksize-8)/4;
  for(unsigned int i = 0; i < numpart; i++){
    ifile.read((char*)&array[i], 4);
  }
  ifile.seekg(4, ios_base::cur);
}

int main(int argc, char** argv){
  // we read a binary file
  ifstream ifile("snapshot_0015", ios::in | ios::binary);
  
  unsigned int numgaspart;
  unsigned int ndark;
  unsigned int nstar;
  char a[5];
  // terminating the char array with a null character allows us to construct
  // strings from it. E.g.
  //   string name(a);
  a[4] = '\0';
  unsigned int blocksize;
  bool flags[4] = {false, false, false, false};
  bool all = false;
  
  // skip header
  // 20: 4 bytes + 4 byte label + 4 byte blocksize + 4 bytes + 4 bytes
  ifile.seekg(20, ios_base::cur);
  ifile.read((char*)&numgaspart, 4);
  ifile.read((char*)&ndark, 4);
  // 8: 2 empty size fields (particle type 2 and 3)
  ifile.seekg(8, ios_base::cur);
  ifile.read((char*)&nstar, 4);
  // 244: the entire header has a size of 256 bytes
  //      20 have gone already
  //      we also have 4 bytes at the end of the block
  //      and 4 bytes at the beginning of the next block
  ifile.seekg(244, ios_base::cur);
  
  float* positions = new float[3*(numgaspart+ndark+nstar)];
  float* densities = new float[numgaspart];
  float* smoothings = new float[numgaspart];
  float* temperatures = new float[numgaspart];
  
  // the !all trick allows us to quit reading early if all arrays have been
  // filled
  while(!all && ifile.good()){
    // read the block label and make a string out of it
    ifile.read(a, 4);
    string name(a);
    bool found = false;
    // check if the string is one of the four labels we need
    if(name == "POS "){
      flags[0] = true;
      found = true;
      read_positions(positions, ifile);
    }
    if(name == "RHO "){
      flags[1] = true;
      found = true;
      read_array(densities, ifile);
    }
    if(name == "HSML"){
      flags[2] = true;
      found = true;
      read_array(smoothings, ifile);
    }
    if(name == "TEMP"){
      flags[3] = true;
      found = true;
      read_array(temperatures, ifile);
    }
    
    // if something was found, we already read the block. If not, we still need
    // to skip it
    if(found){
      // skip beginning of next name block (4 bytes)
      ifile.seekg(4, ios_base::cur);
    } else {
      ifile.read((char*)&blocksize, 4);
      // skip end of name block (4 bytes) + complete block + begin of next name
      // block (4 bytes)
      ifile.seekg(blocksize+8, ios_base::cur);
    }
    all = flags[0] & flags[1] & flags[2] & flags[3];
  }
  
  // do something with the data: output the positions
  for(unsigned int i = 0; i < numgaspart; i++){
    cout << positions[3*i] << "\t" << positions[3*i+1]
         << "\t" << positions[3*i+2] << endl;
  }
  
  delete [] positions;
  delete [] densities;
  delete [] smoothings;
  delete [] temperatures;

  return 0;
}
