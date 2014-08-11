#include <GL/glew.h>
#include <GL/freeglut.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <png.h>
using namespace std;

bool stereo = false;

void parse_options(int argc, char** argv){
  static struct option long_options[] = {
    {"stereo", no_argument, 0, 's'},
    {0, 0, 0, 0}
  };
  
  int c;
  opterr = 0;
  while((c = getopt_long(argc, argv, "s", long_options, NULL)) != -1){
    switch(c){
      case 's':
        stereo = true;
        cout << "Displaying stereo view" << endl;
        break;
      case '?':
        cout << "unknown option: " << argv[optind-1] << endl;
        break;
      default:
        cout << "error" << endl;
        break;
    }
  }
}

string B_to_string(unsigned int bytes){
  stringstream result;
  if(bytes < 1024){
    result << (bytes) << " bytes";
  } else {
    bytes >>= 10;
    if(bytes < 1024){
      result << (bytes) << "KB";
    } else {
      bytes >>= 10;
      result << (bytes) << "MB";
    }
  }
  return result.str();
}

string get_unique_imagename(){
  stringstream fname;
  fname << "image_";
  time_t timev;
  time(&timev);
  fname << timev;
  fname << ".png";
  return fname.str();
}

void unit_matrix(float* matrix){
  for(unsigned int i = 0; i < 16; i++){
    matrix[i] = 0.;
  }
  matrix[0] = 1.;
  matrix[5] = 1.;
  matrix[10] = 1.;
  matrix[15] = 1.;
}

void transpose(float* matrix){
  float newmatrix[16];
  for(unsigned int i = 0; i < 4; i++){
    for(unsigned int j = 0; j < 4; j++){
      newmatrix[4*i+j] = matrix[4*j+i];
    }
  }
  for(unsigned int i = 0; i < 16; i++){
    matrix[i] = newmatrix[i];
  }
}

void multiply(float* A, float* B){
  float C[16] = {0.};
  unsigned int k;
  for(unsigned int i = 0; i < 4; i++){
    for(unsigned int j = 0; j < 4; j++){
      k = 4*i+j;
      for(unsigned int l = 0; l < 4; l++){
        C[k] += A[4*i+l]*B[4*l+j];
      }
    }
  }
  for(unsigned int i = 0; i < 16; i++){
    A[i] = C[i];
  }
}

void premultiply(float* A, float* B){
  float C[16] = {0.};
  unsigned int k;
  for(unsigned int i = 0; i < 4; i++){
    for(unsigned int j = 0; j < 4; j++){
      k = 4*i+j;
      for(unsigned int l = 0; l < 4; l++){
        C[k] += A[4*i+l]*B[4*l+j];
      }
    }
  }
  for(unsigned int i = 0; i < 16; i++){
    B[i] = C[i];
  }
}

void rotate_x(float* matrix, double angle){
  float xrot[16] = {cos(angle), 0., -sin(angle), 0.,
                    0., 1., 0., 0.,
                    sin(angle), 0., cos(angle), 0.,
                    0., 0., 0., 1.};
  multiply(matrix, xrot);
}

void rotate_y(float* matrix, double angle){
  float yrot[16] = {1., 0., 0., 0.,
                    0., cos(angle), sin(angle), 0.,
                    0., -sin(angle), cos(angle), 0.,
                    0., 0., 0., 1.};
  multiply(matrix, yrot);
}

class GLProgram{
private:
  GLuint _program;
  map<string, GLint> _uniforms;
  map<string, GLint> _attributes;

  GLuint get_shader(const char* filename, GLenum shadertype){
    GLuint shader = glCreateShader(shadertype);
    string line, text;
    ifstream in(filename);
    while(getline(in, line)){
      text += line + "\n";
    }
    const char* filecontent = text.c_str();
    glShaderSource(shader, 1, &filecontent, NULL);
    glCompileShader(shader);
    GLint result;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &result);
    if(!result){
      cerr << "Error while compiling shader \"" << filename << "\":\n";
      GLchar infolog[1000];
      GLsizei loglength;
      glGetShaderInfoLog(shader, 1000, &loglength, infolog);
      cerr << infolog << endl;
      exit(1);
    }
    return shader;
  }

  GLuint get_program(const char* vertname, const char* fragname){
    GLuint program = glCreateProgram();
    GLuint vertex_shader = get_shader(vertname, GL_VERTEX_SHADER);
    GLuint fragment_shader = get_shader(fragname, GL_FRAGMENT_SHADER);
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);
    GLint result;
    glGetProgramiv(program, GL_LINK_STATUS, &result);
    if(!result){
      cerr << "Error while linking GPU program:\n";
      GLchar infolog[1000];
      GLsizei loglength;
      glGetProgramInfoLog(program, 1000, &loglength, infolog);
      cerr << infolog;
      exit(1);
    }
    return program;
  }

public:
  GLProgram(const char* vertname, const char* fragname, vector<string> uniforms, vector<string> attributes){
    _program = get_program(vertname, fragname);
    for(unsigned int i = 0; i < uniforms.size(); i++){
      GLint pointer = glGetUniformLocation(_program, uniforms[i].c_str());
      if(pointer < 0){
        cerr << "Requested unknown uniform: \"" << uniforms[i] << "\"" << endl;
        exit(1);
      }
      _uniforms[uniforms[i]] = pointer;
    }
    for(unsigned int i = 0; i < attributes.size(); i++){
      GLint pointer = glGetAttribLocation(_program, attributes[i].c_str());
      if(pointer < 0){
        cerr << "Requested unknown attribute: \"" << attributes[i] << "\"" << endl;
        exit(1);
      }
      glEnableVertexAttribArray(pointer);
      _attributes[attributes[i]] = pointer;
    }
  }
  
  GLint get_attribute(string name){
    return _attributes[name];
  }
  
  GLint get_uniform(string name){
    return _uniforms[name];
  }
  
  void set_uniform(string name, float value){
    glUniform1f(_uniforms[name], value);
  }
  
  GLuint get_program(){
    return _program;
  }
};

class GLComponent{
public:
  virtual void draw(GLProgram* program, float* camera)=0;
  virtual void set_rebuffer()=0;
  virtual void change_strength(float amount)=0;
};

void* window_instance = NULL;

class Window{
private:
  int _dimensions[2];
  
  bool _mousepressed;
  int _mousepos[2];
  
  float _center[3];
  float _camera[3];
  float _rmatrix[16];
  float _cmatrix[16];
  
  GLuint _fb;
  GLuint _rb;
  
  GLProgram* _program;
  vector<GLComponent*> _components;
  
  void rotate(double xangle, double yangle){
    double c1 = cos(xangle);
    double s1 = sin(xangle);
    double c2 = cos(yangle);
    double s2 = sin(yangle);
    float newrmatrix[16] = {c1, 0., s1, 0.,
                            s1*s2, c2, -c1*s2, 0.,
                            -s1*c2, s2, c1*c2, 0.,
                            0., 0., 0., 1.};
    multiply(_rmatrix, newrmatrix);
    
    // found after heavy trial and error; seems to work
    float newcmatrix[16] = {c1, 0., -s1, 0.,
                            -s1*s2, c2, -c1*s2, 0.,
                            s1*c2, s2, c1*c2, 0.,
                            0., 0., 0., 1.};
    transpose(newcmatrix);
    premultiply(newcmatrix, _cmatrix);
    float newcamera[3] = {0., 0., 5.};
    for(unsigned int i = 0; i < 3; i++){
      _camera[i] = 0.;
      for(unsigned int j = 0; j < 3; j++){
        _camera[i] += _cmatrix[4*j+i]*newcamera[j];
      }
    }
    for(unsigned int i = 0; i < _components.size(); i++){
      _components[i]->set_rebuffer();
    }
    glutPostRedisplay();
  }
  
  // since the angle for stereo viewing is always small, we do not recalculate the order for every eye
  // this makes a huge difference in computation time
  void rotate_x(double xangle){
    double c1 = cos(xangle);
    double s1 = sin(xangle);
    float newrmatrix[16] = {c1, 0., s1, 0.,
                            0., 1., 0., 0.,
                            -s1, 0., c1, 0.,
                            0., 0., 0., 1.};
    multiply(_rmatrix, newrmatrix);
    glutPostRedisplay();
  }
  
  void setup_framebuffer(unsigned int* dimensions){
    glGenFramebuffersEXT(1, &_fb);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _fb);
    glGenRenderbuffersEXT(1, &_rb);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _rb);
    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_RGBA8, dimensions[0], dimensions[1]);
    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_RENDERBUFFER_EXT, _rb);
    GLenum status;
    status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    if(status != GL_FRAMEBUFFER_COMPLETE_EXT){
      cerr << "Error in creating image framebuffer" << endl;
    }
    // set up the view with the correct framebuffer dimensions
    glViewport(0, 0, dimensions[0], dimensions[1]);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45., float(dimensions[0])/float(dimensions[1]), 0.1, 100.);
  }
  
  void unset_framebuffer(){
    // reset the view
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    glViewport(0, 0, _dimensions[0], _dimensions[1]);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45., float(_dimensions[0])/float(_dimensions[1]), 0.1, 100.);
    // delete the framebuffer
    glDeleteRenderbuffersEXT(1, &_rb);
    glDeleteFramebuffersEXT(1, &_fb);
  }

public:
  Window(int* argc, char** argv, int width, int height){
    window_instance = this;
    
    _dimensions[0] = width;
    _dimensions[1] = height;
    
    _mousepressed = false;
    _mousepos[0] = 0;
    _mousepos[1] = 0;
    
    _center[0] = 0.;
    _center[1] = 0.;
    _center[2] = -5.;
    
    _camera[0] = 0.;
    _camera[1] = 0.;
    _camera[2] = 5.;
    unit_matrix(_rmatrix);
    unit_matrix(_cmatrix);
    
    glutInit(argc, argv);
    if(stereo){
      glutInitDisplayMode(GLUT_STEREO | GL_DOUBLE);
    } else {
      glutInitDisplayMode(GL_DOUBLE);
    }
    
    glutInitWindowSize(width, height);
    glutCreateWindow("3D SPH viewer");
    glewInit();
    glClearColor(0.,0.,0.,1.);
    
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45., 1., 0.1, 100.);
    
    glutDisplayFunc(displayWrapper);
    glutReshapeFunc(resizeWrapper);
    glutMouseFunc(mousepressWrapper);
    glutMotionFunc(mousemoveWrapper);
    glutSpecialFunc(arrowpressWrapper);
    glutKeyboardFunc(keyfunctionsWrapper);
  }
  
  void add_component(GLComponent* component){
    _components.push_back(component);
  }
  
  void start(GLProgram* program){
    _program = program;
    glutMainLoop();
  }
  
  void displaywrapper(){
    if(stereo){
      glDrawBuffer(GL_BACK_LEFT);
      display();
      double view_angle = 0.034906585;
      rotate_x(view_angle);
      glDrawBuffer(GL_BACK_RIGHT);
      display();
      rotate_x(-view_angle);
    } else {
      display();
    }
    glutSwapBuffers();
  }
  
  void display(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glUseProgram(_program->get_program());
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(_center[0], _center[1], _center[2]);
    glMultMatrixf(_rmatrix);
    for(unsigned int i = 0; i < _components.size(); i++){
      _components[i]->draw(_program, _camera);
    }
  }
  
  void resize(int width, int height){
    _dimensions[0] = width;
    _dimensions[1] = height;
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45., float(width)/float(height), 0.1, 100.);
  }
  
  void mousepress(int button, int state, int x, int y){
    if(button == GLUT_LEFT_BUTTON){
      if(state == GLUT_DOWN){
        _mousepressed = true;
        _mousepos[0] = x;
        _mousepos[1] = y;
      } else {
        _mousepressed = false;
      }
    }
    if(button > 2){
      if(button == 3){
        _center[2] += 0.1;
      } else {
        _center[2] -= 0.1;
      }
      glutPostRedisplay();
    }
  }
  
  void mousemove(int x, int y){
    if(_mousepressed){
      double deltax = M_PI*(x - _mousepos[0])/1800.;
      double deltay = M_PI*(y - _mousepos[1])/1800.;
      _mousepos[0] = x;
      _mousepos[1] = y;
      rotate(-deltax, -deltay);
    }
  }
  
  void arrowpress(int key, int x, int y){
    if(key == GLUT_KEY_UP || key == GLUT_KEY_DOWN){
      if(key == GLUT_KEY_UP){
        _center[1] -= 0.1;
      } else {
        _center[1] += 0.1;
      }
      glutPostRedisplay();
    }
    if(key == GLUT_KEY_LEFT || key == GLUT_KEY_RIGHT){
      if(key == GLUT_KEY_LEFT){
        _center[0] += 0.1;
      } else {
        _center[0] -= 0.1;
      }
      glutPostRedisplay();
    }
  }
  
  void saveImage(unsigned int* dimensions){
    display();
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    char* data = new char[3*dimensions[0]*dimensions[1]];
    glReadPixels(0, 0, dimensions[0], dimensions[1],
                 GL_RGB, GL_UNSIGNED_BYTE, data);
    
    string fname = get_unique_imagename();
    FILE* file = fopen(fname.c_str(), "wb");
    if(!file){
      cerr << "Error: file could not be created!" << endl;
      exit(1);
    }
    
    png_structp png_ptr = 
      png_create_write_struct(PNG_LIBPNG_VER_STRING,
                              (png_voidp)NULL,
                              NULL, NULL);
    if(!png_ptr){
      cerr << "Error allocating png write struct!" << endl;
      exit(1);
    }
  
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if(!info_ptr){
      png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
      cerr << "Error allocating info struct!" << endl;
      exit(1);
    }
  
    if(setjmp(png_jmpbuf(png_ptr))){
      png_destroy_write_struct(&png_ptr, &info_ptr);
      fclose(file);
      cerr << "Error writing PNG" << endl;
      exit(1);
    }
  
    png_init_io(png_ptr, file);
    png_set_IHDR(png_ptr, info_ptr, dimensions[0], dimensions[1], 8,
                 PNG_COLOR_TYPE_RGB,
                 PNG_INTERLACE_ADAM7, PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);
 
    char** rows = new char*[dimensions[1]];
    for(unsigned int i = 0; i < dimensions[1]; i++){
      // luckily, glReadPixels also returns data row per row
      // unfortunately, the offset is the lower left corner, while we start
      // in the upper left corner for our file
      rows[i] = &data[(dimensions[1] - i - 1)*dimensions[0]*3];
    }
    png_set_rows(png_ptr, info_ptr, (png_bytepp)rows);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    png_write_end(png_ptr, info_ptr);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(file);
  
    delete [] rows;
    delete [] data;
    
    cout << "Wrote image file \"" << fname << "\"" << endl;
  }
  
  void change_strength(float amount){
    _components[1]->change_strength(amount);
    glutPostRedisplay();
  }
  
  void keyfunctions(unsigned char key, int x, int y){
    // image dimensions
    unsigned int dimensions[2] = {_dimensions[0]*2, _dimensions[1]*2};
    switch(key){
      case 'p':
        // initialize a separate framebuffer to render to
        setup_framebuffer(dimensions);
        saveImage(dimensions);
        unset_framebuffer();
        break;
      case 'e':
        change_strength(+0.1);
        break;
      case 'd':
        change_strength(-0.1);
        break;
    }
  }
  
  static void displayWrapper(){
    ((Window*)window_instance)->displaywrapper();
  }
  
  static void resizeWrapper(int width, int height){
    ((Window*)window_instance)->resize(width, height);
  }
  
  static void mousepressWrapper(int button, int state, int x, int y){
    ((Window*)window_instance)->mousepress(button, state, x, y);
  }
  
  static void mousemoveWrapper(int x, int y){
    ((Window*)window_instance)->mousemove(x, y);
  }
  
  static void arrowpressWrapper(int key, int x, int y){
    ((Window*)window_instance)->arrowpress(key, x, y);
  }
  
  static void keyfunctionsWrapper(unsigned char key, int x, int y){
    ((Window*)window_instance)->keyfunctions(key, x, y);
  }
};

class Kernel{
private:
  double _h;
  double _norm;
public:
  Kernel(double h, double norm){
    _h = h;
    _norm = norm;
  }
  
  GLuint get_texture(){
    unsigned char* pixels = new unsigned char[256*256];
    for(unsigned int i = 0; i < 256; i++){
      for(unsigned int j = 0; j < 256; j++){
        double x = i*0.00390625 - 0.5;
        double y = j*0.00390625 - 0.5;
        unsigned int offset = i*256+j;
        double r = sqrt(x*x+y*y);
        double u = r/_h;
        if(u > 1.){
          pixels[offset] = 0;
        } else {
          double kernel;
          if(u < 0.5){
            double u2 = u*u;
            kernel = 1. - 6.*u2 + 6.*u*u2;
          } else {
            u = 1.-u;
            u *= u*u;
            kernel = 2.*u;
          }
          pixels[offset] = _norm*kernel*255;
        }
      }
    }
    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 256, 256, 0, GL_RED, GL_UNSIGNED_BYTE, pixels);
    delete [] pixels;
    return texture;
  }
  
  unsigned int get_texture_size(){
    return 256*256*sizeof(unsigned char);
  }
};

class PointSorter{
private:
  double* _depths;

public:
  PointSorter(double* depths){
    _depths = depths;
  }
  
  bool operator()(unsigned int i, unsigned int j){
    return _depths[i] > _depths[j];
  }
};

class PointSet : public GLComponent{
private:
  unsigned int _numpoints;
  unsigned int _numbuffer;
  unsigned int _piecesize;
  unsigned int* _numvertices;
  unsigned int _memsize;
  
  double* _positions;
  double* _smoothings;
  double* _densities;
  double* _temperatures;
  
  float _strength;
  float _density_strength;
  float _max_temperature;
  
  GLuint _texture;
  
  GLuint* _vbo;
  GLuint* _tbo;
  GLuint* _sbo;
  GLuint* _dbo;
  GLuint* _hbo;
  GLuint* _ibo;
  
  bool _do_rebuffer;
  
  void rebuffer(float* camera){
    unsigned int* sorted = new unsigned int[_numpoints];
    double b2 = 0.;
    for(unsigned int j = 0; j < 3; j++){
      b2 += camera[j]*camera[j];
    }
    double* depths = new double[_numpoints];
    for(unsigned int i = 0; i < _numpoints; i++){
      sorted[i] = i;
      double bc = 0.;
      for(unsigned int j = 0; j < 3; j++){
        bc += camera[j]*_positions[3*i+j];
      }
      depths[i] = b2-bc;
    }
    PointSorter sorter(depths);
    sort(&sorted[0], &sorted[_numpoints], sorter);
    for(unsigned int nbuf = 0; nbuf < _numbuffer; nbuf++){
      unsigned int num_this_buffer = min(_piecesize, _numpoints - nbuf*_piecesize);
      float* vertices = new float[12*num_this_buffer];
      float* texcoords = new float[8*num_this_buffer];
      float* allsmoothings = new float[4*num_this_buffer];
      float* alldensities = new float[4*num_this_buffer];
      float* alltemperatures = new float[4*num_this_buffer];
      unsigned short* ids = new unsigned short[4*num_this_buffer];
      _memsize += 36*num_this_buffer*sizeof(float) + 4*num_this_buffer*sizeof(unsigned short);
      // fill buffers
      for(unsigned int i = 0; i < num_this_buffer; i++){
        // we process the points with the greatest depth first
        unsigned int index = sorted[nbuf*_piecesize+i];
        for(unsigned int j = 0; j < 4; j++){
          vertices[12*i+3*j] = _positions[index*3];
          vertices[12*i+3*j+1] = _positions[index*3+1];
          vertices[12*i+3*j+2] = _positions[index*3+2];
          allsmoothings[4*i+j] = _smoothings[index];
          alldensities[4*i+j] = _densities[index];
          alltemperatures[4*i+j] = _temperatures[index];
          ids[4*i+j] = 4*i+j;
        }
        texcoords[8*i] = 0.;
        texcoords[8*i+1] = 0.;
        texcoords[8*i+2] = 1.;
        texcoords[8*i+3] = 0.;
        texcoords[8*i+4] = 1.;
        texcoords[8*i+5] = 1.;
        texcoords[8*i+6] = 0.;
        texcoords[8*i+7] = 1.;
      }
      
      // send data to GPU
      glBindBuffer(GL_ARRAY_BUFFER, _vbo[nbuf]);
      glBufferData(GL_ARRAY_BUFFER, 12*num_this_buffer*sizeof(float), vertices, GL_DYNAMIC_DRAW);
      glBindBuffer(GL_ARRAY_BUFFER, _tbo[nbuf]);
      glBufferData(GL_ARRAY_BUFFER, 8*num_this_buffer*sizeof(float), texcoords, GL_DYNAMIC_DRAW);
      glBindBuffer(GL_ARRAY_BUFFER, _sbo[nbuf]);
      glBufferData(GL_ARRAY_BUFFER, 4*num_this_buffer*sizeof(float), allsmoothings, GL_DYNAMIC_DRAW);
      glBindBuffer(GL_ARRAY_BUFFER, _dbo[nbuf]);
      glBufferData(GL_ARRAY_BUFFER, 4*num_this_buffer*sizeof(float), alldensities, GL_DYNAMIC_DRAW);
      glBindBuffer(GL_ARRAY_BUFFER, _hbo[nbuf]);
      glBufferData(GL_ARRAY_BUFFER, 4*num_this_buffer*sizeof(float), alltemperatures, GL_DYNAMIC_DRAW);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _ibo[nbuf]);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, 4*num_this_buffer*sizeof(unsigned short), ids, GL_DYNAMIC_DRAW);
      
      _numvertices[nbuf] = 4*num_this_buffer;
      
      delete [] vertices;
      delete [] texcoords;
      delete [] allsmoothings;
      delete [] alldensities;
      delete [] alltemperatures;
      delete [] ids;
    }
    delete [] sorted;
    delete [] depths;
  }

public:
  PointSet(vector<double>& positions, vector<double>& smoothings, vector<double>& densities, vector<double>& temperatures, Kernel& kernel,
           float strength = 1., float density_strength = 1., float max_temperature = 1.){
    _numpoints = smoothings.size();
    // the second part ceils the integer result
    _numbuffer = 4*_numpoints/65536 + ((4*_numpoints)%65536 > 0);
    // integer division automatically floors, which is what we want here
    _piecesize = _numpoints/_numbuffer;
    _memsize = 0;
    _numvertices = new unsigned int[_numbuffer];
    
    _positions = new double[3*_numpoints];
    _smoothings = new double[_numpoints];
    _densities = new double[_numpoints];
    _temperatures = new double[_numpoints];
    for(unsigned int i = 0; i < _numpoints; i++){
      for(unsigned int j = 0; j < 3; j++){
        _positions[3*i+j] = positions[3*i+j];
      }
      _smoothings[i] = smoothings[i];
      _densities[i] = densities[i];
      _temperatures[i] = temperatures[i];
    }
    
    _strength = strength;
    _density_strength = density_strength;
    _max_temperature = max_temperature;
    
    _texture = kernel.get_texture();
    _memsize += kernel.get_texture_size();
    
    _vbo = new GLuint[_numbuffer];
    _tbo = new GLuint[_numbuffer];
    _sbo = new GLuint[_numbuffer];
    _dbo = new GLuint[_numbuffer];
    _hbo = new GLuint[_numbuffer];
    _ibo = new GLuint[_numbuffer];
    glGenBuffers(_numbuffer, _vbo);
    glGenBuffers(_numbuffer, _tbo);
    glGenBuffers(_numbuffer, _sbo);
    glGenBuffers(_numbuffer, _dbo);
    glGenBuffers(_numbuffer, _hbo);
    glGenBuffers(_numbuffer, _ibo);
    float camera[3] = {0., 0., 5.};
    rebuffer(camera);
    _do_rebuffer = false;
    cout << "Allocated " << B_to_string(_memsize) << " on GPU" << endl;
  }
  
  ~PointSet(){
    delete [] _positions;
    delete [] _smoothings;
    delete [] _densities;
    delete [] _temperatures;
    
    delete [] _numvertices;
    delete [] _vbo;
    delete [] _tbo;
    delete [] _sbo;
    delete [] _dbo;
    delete [] _hbo;
    delete [] _ibo;
  }
  
  void draw(GLProgram* program, float* camera){
    if(_do_rebuffer){
      rebuffer(camera);
      _do_rebuffer = false;
    }
    program->set_uniform("strength", _strength);
    program->set_uniform("densityStrength", _density_strength);
    program->set_uniform("maxTemperature", _max_temperature);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, _texture);
    glUniform1i(program->get_uniform("uSampler"), 0);
    for(unsigned int nbuf = 0; nbuf < _numbuffer; nbuf++){
      glBindBuffer(GL_ARRAY_BUFFER, _vbo[nbuf]);
      glVertexAttribPointer(program->get_attribute("aVertexPosition"), 3, GL_FLOAT, GL_FALSE, 0, NULL);
      glBindBuffer(GL_ARRAY_BUFFER, _tbo[nbuf]);
      glVertexAttribPointer(program->get_attribute("aTextureCoord"), 2, GL_FLOAT, GL_FALSE, 0, NULL);
      glBindBuffer(GL_ARRAY_BUFFER, _sbo[nbuf]);
      glVertexAttribPointer(program->get_attribute("aSmoothing"), 1, GL_FLOAT, GL_FALSE, 0, NULL);
      glBindBuffer(GL_ARRAY_BUFFER, _dbo[nbuf]);
      glVertexAttribPointer(program->get_attribute("aDensity"), 1, GL_FLOAT, GL_FALSE, 0, NULL);
      glBindBuffer(GL_ARRAY_BUFFER, _hbo[nbuf]);
      glVertexAttribPointer(program->get_attribute("aTemperature"), 1, GL_FLOAT, GL_FALSE, 0, NULL);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _ibo[nbuf]);
      glDrawElements(GL_QUADS, _numvertices[nbuf], GL_UNSIGNED_SHORT, NULL);
    }
  }
  
  void set_rebuffer(){
    _do_rebuffer = true;
  }
  
  void change_strength(float amount){
    _density_strength += amount;
  }
};

void read_positions(float* positions, ifstream& ifile){
  unsigned int blocksize;
  ifile.read((char*)&blocksize, 4);
  // skip end of name block and beginning of data block
  ifile.seekg(8, ios_base::cur);
  unsigned int numpart = (blocksize-8)/12;
  for(unsigned int i = 0; i < numpart; i++){
    ifile.read((char*)&positions[3*i], 12);
  }
  ifile.seekg(4, ios_base::cur);
}

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
  parse_options(argc, argv);
  Window window(&argc, argv, 400, 400);
  
  vector<string> uniforms;
  uniforms.push_back("strength");
  uniforms.push_back("uSampler");
  uniforms.push_back("densityStrength");
  uniforms.push_back("maxTemperature");
  vector<string> attributes;
  attributes.push_back("aVertexPosition");
  attributes.push_back("aTextureCoord");
  attributes.push_back("aSmoothing");
  attributes.push_back("aDensity");
  attributes.push_back("aTemperature");
  GLProgram program("shader.vert", "shader.frag", uniforms, attributes);
  
  ifstream ifile("snapshot_0015", ios::in | ios::binary);
  
  unsigned int numgaspart;
  unsigned int ndark;
  unsigned int nstar;
  char a[5];
  a[4] = '\0';
  unsigned int blocksize;
  bool flags[4] = {false, false, false, false};
  bool all = false;
  
  // skip header
  ifile.seekg(20, ios_base::cur);
  ifile.read((char*)&numgaspart, 4);
  ifile.read((char*)&ndark, 4);
  ifile.seekg(8, ios_base::cur);
  ifile.read((char*)&nstar, 4);
  ifile.seekg(244, ios_base::cur);
  
  float* positions = new float[3*(numgaspart+ndark+nstar)];
  float* densities = new float[numgaspart];
  float* smoothings = new float[numgaspart];
  float* temperatures = new float[numgaspart];
  
  while(!all && ifile.good()){
    ifile.read(a, 4);
    string name(a);
    bool found = false;
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
  
  vector<double> gaspositions(3*numgaspart);
  vector<double> gassmoothings(numgaspart);
  vector<double> gasdensities(numgaspart);
  vector<double> gastemperatures(numgaspart);
  double maxgasdensity = 0.;
  double maxgastemperature = 2.e4;
  for(unsigned int i = 0; i < numgaspart; i++){
    gaspositions[3*i] = positions[3*i];
    gaspositions[3*i+1] = positions[3*i+1];
    gaspositions[3*i+2] = positions[3*i+2];
    gasdensities[i] = densities[i];
    gassmoothings[i] = smoothings[i];
    gastemperatures[i] = temperatures[i];
    maxgasdensity = max(maxgasdensity, gasdensities[i]);
  }
  for(unsigned int i = 0; i < numgaspart; i++){
    gasdensities[i] /= maxgasdensity;
    gastemperatures[i] /= maxgastemperature;
  }
  Kernel kernel(0.5, 32./3./50.);
  PointSet gas(gaspositions, gassmoothings, gasdensities, gastemperatures, kernel);
  
  vector<double> starpositions(3*nstar);
  vector<double> starsmoothings(nstar);
  vector<double> stardensities(nstar);
  vector<double> startemperatures(nstar);
  unsigned int offset = 3*(numgaspart+ndark);
  for(unsigned int i = 0; i < nstar; i++){
    starpositions[3*i] = positions[offset+3*i];
    starpositions[3*i+1] = positions[offset+3*i+1];
    starpositions[3*i+2] = positions[offset+3*i+2];
    starsmoothings[i] = 0.025;
    stardensities[i] = 0.00116;
    startemperatures[i] = 0.00116;
  }
  Kernel kernel2(0.04, 0.5);
  PointSet stars(starpositions, starsmoothings, stardensities, startemperatures, kernel2, 100., 100.);
  delete [] positions;
  delete [] densities;
  delete [] smoothings;
  delete [] temperatures;
  
//  ifstream file("gas.dat", ios::in|ios::binary|ios::ate);
//  streampos size = file.tellg();
//  float * data = new float[size/sizeof(float)];
//  file.seekg(0, ios::beg);
//  file.read(reinterpret_cast<char*>(data), size);
//  unsigned int numpoints = size/(6*sizeof(float));
//  vector<double> positions(3*numpoints);
//  vector<double> smoothings(numpoints);
//  vector<double> densities(numpoints);
//  vector<double> temperatures(numpoints);
//  double maxdensity = 0.;
//  double maxtemperature = 2.e4;
//  for(unsigned int i = 0; i < numpoints; i++){
//    positions[3*i] = data[6*i];
//    positions[3*i+1] = data[6*i+1];
//    positions[3*i+2] = data[6*i+2];
//    smoothings[i] = data[6*i+4];
//    densities[i] = data[6*i+3];
//    temperatures[i] = data[6*i+5];
//    maxdensity = max(maxdensity, densities[i]);
//  }
//  delete [] data;
//  file.close();
//  for(unsigned int i = 0; i < numpoints; i++){
//    densities[i] /= maxdensity;
//    temperatures[i] /= maxtemperature;
//  }
//  Kernel kernel(0.5, 32./3./50.);
//  PointSet pointset(positions, smoothings, densities, temperatures, kernel);
//  
//  file.open("stars.dat", ios::in|ios::binary|ios::ate);
//  size = file.tellg();
//  data = new float[size/sizeof(float)];
//  file.seekg(0, ios::beg);
//  file.read(reinterpret_cast<char*>(data), size);
//  numpoints = size/(3*sizeof(float));
//  positions.resize(3*numpoints);
//  smoothings.resize(numpoints);
//  densities.resize(numpoints);
//  temperatures.resize(numpoints);
//  for(unsigned int i = 0; i < numpoints; i++){
//    positions[3*i] = data[3*i];
//    positions[3*i+1] = data[3*i+1];
//    positions[3*i+2] = data[3*i+2];
//    smoothings[i] = 0.025;
//    densities[i] = 0.00116;
//    temperatures[i] = 0.00116;
//  }
//  delete [] data;
//  file.close();
//  Kernel kernel2(0.04, 0.5);
//  PointSet stars(positions, smoothings, densities, temperatures, kernel2, 100., 100.);
  
  window.add_component(&stars);
  window.add_component(&gas);
  
  window.start(&program);
  
  return 0;
}
