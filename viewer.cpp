#include <GL/glew.h>
#include <GL/freeglut.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

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
    glutInitDisplayMode(GL_DOUBLE);
    
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
    glutSwapBuffers();
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
  
  void keyfunctions(unsigned char key, int x, int y){
    
  }
  
  static void displayWrapper(){
    ((Window*)window_instance)->display();
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
};

int main(int argc, char** argv){
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
  
  ifstream file("gas.dat", ios::in|ios::binary|ios::ate);
  streampos size = file.tellg();
  float * data = new float[size/sizeof(float)];
  file.seekg(0, ios::beg);
  file.read(reinterpret_cast<char*>(data), size);
  unsigned int numpoints = size/(6*sizeof(float));
  vector<double> positions(3*numpoints);
  vector<double> smoothings(numpoints);
  vector<double> densities(numpoints);
  vector<double> temperatures(numpoints);
  double maxdensity = 0.;
  double maxtemperature = 2.e4;
  for(unsigned int i = 0; i < numpoints; i++){
    positions[3*i] = data[6*i];
    positions[3*i+1] = data[6*i+1];
    positions[3*i+2] = data[6*i+2];
    smoothings[i] = data[6*i+4];
    densities[i] = data[6*i+3];
    temperatures[i] = data[6*i+5];
    maxdensity = max(maxdensity, densities[i]);
  }
  delete [] data;
  file.close();
  for(unsigned int i = 0; i < numpoints; i++){
    densities[i] /= maxdensity;
    temperatures[i] /= maxtemperature;
  }
  Kernel kernel(0.5, 32./3./50.);
  PointSet pointset(positions, smoothings, densities, temperatures, kernel);
  
  file.open("stars.dat", ios::in|ios::binary|ios::ate);
  size = file.tellg();
  data = new float[size/sizeof(float)];
  file.seekg(0, ios::beg);
  file.read(reinterpret_cast<char*>(data), size);
  numpoints = size/(3*sizeof(float));
  positions.resize(3*numpoints);
  smoothings.resize(numpoints);
  densities.resize(numpoints);
  temperatures.resize(numpoints);
  for(unsigned int i = 0; i < numpoints; i++){
    positions[3*i] = data[3*i];
    positions[3*i+1] = data[3*i+1];
    positions[3*i+2] = data[3*i+2];
    smoothings[i] = 0.025;
    densities[i] = 0.00116;
    temperatures[i] = 0.00116;
  }
  delete [] data;
  file.close();
  Kernel kernel2(0.04, 0.5);
  PointSet stars(positions, smoothings, densities, temperatures, kernel2, 100., 100.);
  
  window.add_component(&stars);
  window.add_component(&pointset);
  
  window.start(&program);
  
  return 0;
}
