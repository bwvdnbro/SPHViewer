/*******************************************************************************
 * This file is part of SPHViewer
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * SPHViewer is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SPHViewer is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with SPHViewer. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#include <iostream>
#include <png.h>
#include <cstdlib>
using namespace std;

// Since libpng is not the easiest library to work with: a simple example.
// The following code creates an image file test.png which contains a
// vertical gradient from almost completely transparent to white.
// The image is 200x200 pixels in size and has a 16-bit color depth in RGBA
//
// To compile this program, don't forget to include libpng:
//   g++ -o pngtest pngtest.cpp -lpng

int main(int argc, char** argv){
  // libpng does not do I/O; the file has to be opened separately
  FILE* file = fopen("test.png", "wb");
  if(!file){
    cerr << "Error: file could not be created!" << endl;
    exit(1);
  }
  
  // libpng needs two structs before even starting to work:
  //  - a png_structp which basically identifies the file you will write
  //  - a png_infop which contains information about the image (resolution...)
  // we create the png_structp with empty custom error and warning handle
  // functions
  png_structp png_ptr = 
    png_create_write_struct(PNG_LIBPNG_VER_STRING,
                            (png_voidp)NULL,
                            NULL, NULL);
  if(!png_ptr){
    cerr << "Error allocating png write struct!" << endl;
    exit(1);
  }
  
  // the png_infop requires no specific arguments
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if(!info_ptr){
    png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    cerr << "Error allocating info struct!" << endl;
    exit(1);
  }
  
  // function that is called if an error occurs during the file write
  // we clean up and display an error message
  if(setjmp(png_jmpbuf(png_ptr))){
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(file);
    cerr << "Error writing PNG" << endl;
    exit(1);
  }
  
  // this creates an empty file
  png_init_io(png_ptr, file);
  // set the image info:
  // width, height, bit depth, color type, interlace type, compression type,
  // filter method
  png_set_IHDR(png_ptr, info_ptr, 200, 200, 16, PNG_COLOR_TYPE_RGB_ALPHA,
               PNG_INTERLACE_ADAM7, PNG_COMPRESSION_TYPE_DEFAULT,
               PNG_FILTER_TYPE_DEFAULT);
  
  // initialize image data
  // we need an array with 200 rows, which in turn contain 2*4*200 chars
  // (200 columns, 4 RGBA color components, 16 bits = 2 bytes per component)
  char** rows = new char*[200];
  for(unsigned int i = 0; i < 200; i++){
    rows[i] = new char[1600];
    for(unsigned int j = 0; j < 1600; j++){
      rows[i][j] = i+55;
    }
  }
  // tell libpng about the rows
  png_set_rows(png_ptr, info_ptr, (png_bytepp)rows);
  
  // do the actual write without applying a transformation to the data
  png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
  
  // we are done writing
  png_write_end(png_ptr, info_ptr);
  
  // it's nice to clean up once you're done
  png_destroy_write_struct(&png_ptr, &info_ptr);
  fclose(file);
  
  for(unsigned int i = 0; i < 200; i++){
    delete [] rows[i];
  }
  delete [] rows;
  
  return 0;
}
