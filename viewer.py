################################################################################
# This file is part of SPHViewer
# Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# SPHViewer is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SPHViewer is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with SPHViewer. If not, see <http://www.gnu.org/licenses/>.
################################################################################

from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
import numpy as np
import sys
import struct
import subprocess
import time

# set this variable to enable stereo rendering (for 3D screens)
show_textbox = True

# range of snapshots to be displayed (requires gas_XXXX.dat and stars_XXXX.dat files to be present in this range)
snaps = [169, 169]

## Read and compile a shader of the given type
#
# The source code for the shader is read in from the given filename
# Supported types are GL_VERTEX_SHADER and GL_FRAGMENT_SHADER
#
# @param filename Name of the file containing the raw GSLS source of the shader
# @param shadertype Type of the shader (either GL_VERTEX_SHADER or GL_FRAGMENT_SHADER)
# @return A GL id pointing to the created shader
def get_shader(filename, shadertype):
    shader = glCreateShader(shadertype)
    glShaderSource(shader, open(filename).read())
    glCompileShader(shader)
    result = glGetShaderiv(shader, GL_COMPILE_STATUS)
    if not(result):
      raise RuntimeError(glGetShaderInfoLog(shader))
    return shader

## Compile the given shaders and link them into a GPU program
#
# The shaders are read in from the files with the given names
#
# @param vertname Name of the file containing the raw GSLS source code for the vertex shader
# @param fragname Name of the file containing the raw GSLS source code for the fragment shader
# @return A GL id pointing to the created program
def get_program(vertname, fragname):
    program = glCreateProgram()
    vertex_shader = get_shader(vertname, GL_VERTEX_SHADER)
    fragment_shader = get_shader(fragname, GL_FRAGMENT_SHADER)
    glAttachShader(program, vertex_shader)
    glAttachShader(program, fragment_shader)
    glLinkProgram(program)
    result = glGetProgramiv(program, GL_LINK_STATUS)
    if not(result):
      raise RuntimeError(glGetProgramInfoLog(program))
    return program

## Wrapper around GPU program
#
# Contains the GL id of the program on the GPU and GL id's pointing to the
# uniform variables and attributes on the GPU
# The variable and attribute names are passed on as a list to the constructor and are
# stored as dictionaries, using the respective names as keys
class Program:
  ## Constructor.
  #
  # Compiles and links the program and sets the pointers to the variables and attributes.
  #
  # @param self The pointer to this object
  # @param vertname Name of the file containing the GSLS source code for the vertex shader
  # @param fragname Name of the file containing the GSLS source code for the fragment shader
  # @param uniforms List containing the names of uniform variables in the shaders that will be accessed
  # @param attributes List containing the names of attributes in the shaders that will be accessed
  def __init__(self, vertname, fragname, uniforms, attributes):
    self._program = get_program(vertname, fragname)
    self._uniforms = {}
    for uniform in uniforms:
      self._uniforms[uniform] = glGetUniformLocation(self._program, uniform)
    self._attributes = {}
    for attribute in attributes:
      self._attributes[attribute] = glGetAttribLocation(self._program, attribute)
      glEnableVertexAttribArray(self._attributes[attribute])
    
    self._pointsets = {}
  
  ## Add a pointset to the program
  #
  # Pointsets are stored in a dictionary for easy access
  #
  # @param self The pointer to this object
  # @param name Key for the pointset in the internal dictionary
  # @param pointset Pointset instance
  def add_pointset(self, name, pointset):
    self._pointsets[name] = pointset

## Representation of a set of points
#
# Upon construction, the apropriate buffers are initialized on the GPU
# The pointset can be drawn by calling its draw function
class PointSet:
  ## Constructor.
  #
  # The data are split up into GPU manageable parts (4-byte indices = arrays of at most 65536 elements)
  # For every chunck of data, 6 buffers are constructed, containing
  #  * vertex positions
  #  * texture coordinates
  #  * smoothing lengths
  #  * densities
  #  * temperatures
  #  * element indices to pass on to glDrawElements()
  # Every point is represented as a square and hence corresponds to 4 values in every buffer
  # Vertex positions are 3-dimensional, while texture coordinates are 2-dimensional
  #
  # @param self The pointer to this object
  # @param positions Array of length N containing triples of coordinates
  # @param smoothings Array of length N containing smoothing lengths
  # @param densities Array of length N containing densities
  # @param temperatures Array of length N containing temperatures
  def __init__(self, positions, smoothings, densities, temperatures):
    self._numPoints = len(positions)
    self._numBuffer = max(1,int(np.ceil(4.*self._numPoints/65536.)))
    self._pieceSize = int(np.floor(np.float(self._numPoints)/np.float(self._numBuffer)))
    self._vbo = [0]*self._numBuffer
    self._tbo = [0]*self._numBuffer
    self._sbo = [0]*self._numBuffer
    self._dbo = [0]*self._numBuffer
    self._hbo = [0]*self._numBuffer
    self._ibo = [0]*self._numBuffer
    self._numtriangles = [0]*self._numBuffer
    self._memsize = 0
    for nbuf in range(self._numBuffer):
      numThisBuffer = min(self._pieceSize, self._numPoints - nbuf*self._pieceSize)
      vertices = np.zeros(12*numThisBuffer, "f4")
      texcoords = np.zeros(8*numThisBuffer, "f4")
      allsmoothings = np.zeros(4*numThisBuffer, "f4")
      alldensities = np.zeros(4*numThisBuffer, "f4")
      alltemperatures = np.zeros(4*numThisBuffer, "f4")
      ids = np.zeros(4*numThisBuffer, "u2")
      self._memsize += vertices.nbytes + texcoords.nbytes + allsmoothings.nbytes + alldensities.nbytes + alltemperatures.nbytes + ids.nbytes
      for i in range(numThisBuffer):
        for j in range(4):
          vertices[12*i+3*j] = positions[nbuf*numThisBuffer+i][0]
          vertices[12*i+3*j+1] = positions[nbuf*numThisBuffer+i][1]
          vertices[12*i+3*j+2] = positions[nbuf*numThisBuffer+i][2]
          allsmoothings[4*i+j] = smoothings[nbuf*numThisBuffer+i]
          alldensities[4*i+j] = densities[nbuf*numThisBuffer+i]
          alltemperatures[4*i+j] = temperatures[nbuf*numThisBuffer+i]
          ids[4*i+j] = 4*i+j
        texcoords[8*i] = 0.0
        texcoords[8*i+1] = 0.0
        texcoords[8*i+2] = 1.0
        texcoords[8*i+3] = 0.0
        texcoords[8*i+4] = 1.0
        texcoords[8*i+5] = 1.0
        texcoords[8*i+6] = 0.0
        texcoords[8*i+7] = 1.0
      self._vbo[nbuf] = glGenBuffers(1)
      glBindBuffer(GL_ARRAY_BUFFER, self._vbo[nbuf])
      glBufferData(GL_ARRAY_BUFFER, len(vertices)*4, vertices, GL_STATIC_DRAW)
      self._tbo[nbuf] = glGenBuffers(1)
      glBindBuffer(GL_ARRAY_BUFFER, self._tbo[nbuf])
      glBufferData(GL_ARRAY_BUFFER, len(texcoords)*4, texcoords, GL_STATIC_DRAW)
      self._numtriangles[nbuf] = len(allsmoothings)
      self._sbo[nbuf] = glGenBuffers(1)
      glBindBuffer(GL_ARRAY_BUFFER, self._sbo[nbuf])
      glBufferData(GL_ARRAY_BUFFER, len(allsmoothings)*4, allsmoothings, GL_STATIC_DRAW)
      self._dbo[nbuf] = glGenBuffers(1)
      glBindBuffer(GL_ARRAY_BUFFER, self._dbo[nbuf])
      glBufferData(GL_ARRAY_BUFFER, len(alldensities)*4, alldensities, GL_STATIC_DRAW)
      self._hbo[nbuf] = glGenBuffers(1)
      glBindBuffer(GL_ARRAY_BUFFER, self._hbo[nbuf])
      glBufferData(GL_ARRAY_BUFFER, len(alltemperatures)*4, alltemperatures, GL_STATIC_DRAW)
      self._ibo[nbuf] = glGenBuffers(1)
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self._ibo[nbuf])
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, len(ids)*2, ids, GL_STATIC_DRAW)
  
  ## Draw the vertices in the pointset to the active OpenGL draw buffer
  #
  # We first link all 5 attribute buffers to the corresponding attributes
  # We then link the element buffer and call glDrawElements() on it, which does the drawing
  #
  # @param self Pointer to this object
  # @param glprog Program used to draw on the GPU. Contains pointers to the attribute buffers
  def draw(self, glprog):
    for nbuf in range(self._numBuffer):
      glBindBuffer(GL_ARRAY_BUFFER, self._vbo[nbuf])
      glVertexAttribPointer(glprog._attributes["aVertexPosition"], 3, GL_FLOAT, GL_FALSE, 0, None)
      glBindBuffer(GL_ARRAY_BUFFER, self._tbo[nbuf])
      glVertexAttribPointer(glprog._attributes["aTextureCoord"], 2, GL_FLOAT, GL_FALSE, 0, None)
      glBindBuffer(GL_ARRAY_BUFFER, self._sbo[nbuf])
      glVertexAttribPointer(glprog._attributes["aSmoothing"], 1, GL_FLOAT, GL_FALSE, 0, None)
      glBindBuffer(GL_ARRAY_BUFFER, self._dbo[nbuf])
      glVertexAttribPointer(glprog._attributes["aDensity"], 1, GL_FLOAT, GL_FALSE, 0, None)
      glBindBuffer(GL_ARRAY_BUFFER, self._hbo[nbuf])
      glVertexAttribPointer(glprog._attributes["aTemperature"], 1, GL_FLOAT, GL_FALSE, 0, None)
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self._ibo[nbuf])
      glDrawElements(GL_QUADS, self._numtriangles[nbuf], GL_UNSIGNED_SHORT, None)

## Text line on the bottom left of the screen
#
# The displayed information is stored internally
# To draw the text line, the get_text() method is called
class TextBox:
  ## Constructor.
  #
  # Initialize variables
  #
  # @param self Pointer to this object
  # @param variable Name of the variable that is used for color scaling
  # @param strength Strength with which the opacities are scaled in the fragment shader
  # @param maxTemperature Maximal temperature of the gas (used for the gas color)
  # @param curSnap Index of the currently active snapshot
  def __init__(self, variable, strength, maxTemperature, curSnap):
    self._variable = variable
    self._strength = strength
    self._maxTemperature = maxTemperature
    self._curSnap = curSnap
  
  ## Generate the text line
  #
  # The text line contains formatted information on the current plotted variables
  #
  # @param self Pointer to this object
  # @return String to be displayed in the bottom left corner of the screen
  def get_text(self):
    return self._variable + ", Strength: %.02f" % 10.**self._strength + ", MaxTemperature: %.0f" % self._maxTemperature + ", Current snapshot: %i" % self._curSnap

## Container for information about the current view
#
# Contains the center and rotation of the view and
# the dimensions of the window
class View:
  ## Constructor.
  #
  # Initialize variables
  #
  # @param self The pointer to this object
  # @param center The center of the view ([center_x, center_y, center_z])
  # @param rmatrix Rotation matrix (4x4 numpy array)
  # @param dimensions Dimensions of the window ([width, height])
  # @param stereo Flag signaling if OpenGL stereo is enabled
  def __init__(self, center = [0., 0., -5.], rmatrix = np.identity(4), dimensions = [400, 400], stereo = False):
    self._center = center
    self._rmatrix = rmatrix
    self._dimensions = dimensions
    self._stereo = stereo
  
  ## Translate the center of the view
  #
  # @param self The pointer to this object
  # @param vector Translation vector for the center ([x_translation, y_translation, z_translation])
  def translate(self, vector):
    self._center[0] += vector[0]
    self._center[1] += vector[1]
    self._center[2] += vector[2]
  
  ## Rotate the view
  #
  # The rotation is split in a rotation around the x-axis and one around the y-axis
  #
  # @param self The pointer to this object
  # @param xangle Angle for rotation around the x-axis
  # @param yangle Angle for rotation around the y-axis
  def rotate(self, xangle, yangle):
    newrmatrix = np.identity(4)
    yrot = np.array([[1.,0.,0.,0.],
                     [0., np.cos(yangle/10.), np.sin(yangle/10.), 0.],
                     [0., -np.sin(yangle/10.), np.cos(yangle/10.), 0.],
                     [0., 0., 0., 1.]])
    xrot = np.array([[np.cos(xangle/10.), 0., -np.sin(xangle/10.), 0.],
                     [0., 1., 0., 0.],
                     [np.sin(xangle/10.), 0., np.cos(xangle/10.), 0.],
                     [0., 0., 0., 1.]])
    newrmatrix = np.dot(newrmatrix,xrot)
    newrmatrix = np.dot(newrmatrix,yrot)
    self._rmatrix = np.dot(self._rmatrix, newrmatrix)
  
  def rotate_x(self, xangle):
    newrmatrix = np.identity(4)
    xrot = np.array([[np.cos(xangle), 0., np.sin(xangle), 0.],
                     [0., 1., 0., 0.],
                     [-np.sin(xangle), 0., np.cos(xangle), 0.],
                     [0., 0., 0., 1.]])
    newrmatrix = np.dot(newrmatrix,xrot)
    self._rmatrix = np.dot(self._rmatrix, newrmatrix)
  
  ## Resize the view
  #
  # @param self The pointer to this object
  # @param new_size New dimensions for the view ([width, height])
  def resize(self, new_size):
    self._dimensions = new_size

## global variables, since we have to use a static display function
# GPU program
glprog = None
# texture used to plot gas particles
texture = 0
# texture used to plot stars
startexture = 0
# strength of the scaling used for the opacity in the fragment shader
gasstrength = 0.
# Text line in the bottom left corner of the screen
textbox = None
# Active snapshot
active = snaps[0]
# Maximal temperature
maxTemperature = 4.
# The view. Contains information about the dimensions of the window, the center of the view and its rotation
view = None

## Main function
#
# Constructs a window and a GPU program and starts the main loop for displaying the window
def main():
    global glprog, glname, texture, pointset, stars, startexture, dimensions, textbox, gasstrength, maxTemperature, snaps, view
    
    view = View(stereo = False)
    
    glutInit(sys.argv)
    if view._stereo:
      glutInitDisplayMode(GLUT_STEREO | GL_DOUBLE)
    else:
      glutInitDisplayMode(GL_DOUBLE)
                
    glutInitWindowSize(view._dimensions[0], view._dimensions[1])
    glutCreateWindow("3D SPH viewer")

    glClearColor(0.,0.,0.,1.)
    # about blending:
    #  we use blending to create a transparency effect for the gas: the denser the gas, the less transparent
    #  however, rendering this is order dependant: we can only correctly account for the opacity if we render
    #  the gas particles that are further away first
    #  this requires us to sort the particles. this works fine for a fixed camera position, since then the
    #  sorting has to be done only once. if we want to be able to turn the gas around, this would require
    #  a resort for every rotation. even on a GPU, this is very difficult
    #
    #  more correct would be to use additive blending (setting the second argument below to GL_ONE achieves this)
    #  the problem here is that colors get mixed up: as long as you use one color, additive blending indeed yields
    #  a correct total attenuation of the individual color channels. when you use a combination of color channels,
    #  the resulting color will be a mix of individually attenuated colors. suppose you have a blue gas blob hiding
    #  behind a red gas blob, then both the red and blue gas will be attenuated by the other gas blob, resulting in
    #  slightly opaque red and blue channels. the total color will be white, since red and blue will be equally strong
    #  what you really want in this case is a very vague or invisible blue blob, covered by a strongly visible red blob
    #  but to get this, you have to somehow damp the blue channel, which only works if you know that the blue blob is
    #  behind the red one during rendering
    #
    #  long story short: the blending function below will give some strange behaviour if you rotate the view, but nothing
    #  can be done about this without a lot of extra effort (which would probably slow down the whole application too much)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
    glEnable(GL_BLEND)
    glprog = Program("shader.vert", "shader.frag",
                     ["strength", "uSampler", "densityStrength", "maxTemperature"],
                     ["aVertexPosition", "aTextureCoord", "aSmoothing", "aDensity", "aTemperature"])
    
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(45.0, 1.0, 0.1, 100.0)
    totmemsize = 0
    
    texture = glGenTextures(1)
    glBindTexture(GL_TEXTURE_2D, texture)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE)
    # u1: unsigned 8 bit integers
    pixels = np.zeros(256*256*4, "u1")
    h = 0.5
    norm = 32./3./50.
    for i in range(256):
      for j in range(256):
        x = i*1.0/256.0 - 0.5
        y = j*1.0/256.0 - 0.5
        offset = (i*256+j)*4
        r = np.sqrt(x**2+y**2)
        u = r/h
        if u > 1.0:
          pixels[offset] = 0
          pixels[offset+1] = 0
          pixels[offset+2] = 0
          pixels[offset+3] = 0
        else:
          if u < 0.5:
            u2 = u*u
            kernel = 1.0 - 6.0*u2 + 6.0*u*u2
          else:
            u = 1.0-u
            u *= u*u
            kernel = 2.0*u
          pixels[offset] = 255
          pixels[offset+1] = 0
          pixels[offset+2] = 0
          pixels[offset+3] = norm*kernel*255
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 256, 256, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels)
    totmemsize += pixels.nbytes
    
    startexture = glGenTextures(1)
    glBindTexture(GL_TEXTURE_2D, startexture)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE)
    # u1: unsigned 8 bit integers
    pixels = np.zeros(256*256*4, "u1")
    h = 0.04
    norm = 0.5
    for i in range(256):
      for j in range(256):
        x = i*1.0/256.0 - 0.5
        y = j*1.0/256.0 - 0.5
        offset = (i*256+j)*4
        r = np.sqrt(x**2+y**2)
        u = r/h
        if u > 1.0:
          pixels[offset] = 0
          pixels[offset+1] = 0
          pixels[offset+2] = 0
          pixels[offset+3] = 0
        else:
          if u < 0.5:
            u2 = u*u
            kernel = 1.0 - 6.0*u2 + 6.0*u*u2
          else:
            u = 1.0-u
            u *= u*u
            kernel = 2.0*u
          pixels[offset] = 255
          pixels[offset+1] = 0
          pixels[offset+2] = 0
          pixels[offset+3] = norm*kernel*255
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 256, 256, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels)
    totmemsize += pixels.nbytes
    
    if len(sys.argv) > 1:
      print "generating temporary files"
      command = "hyplot --visual=extract_data.py " + sys.argv[1]
      print command
      hyplot = subprocess.Popen(command, shell=True, executable="/bin/bash", stderr=subprocess.PIPE)
      hyplot.wait()
    
    maxdensity = 0.
    print "calculating maximal density..."
    for snap in range(snaps[0],snaps[1]+1):
      gasparts = open("gas_{0:04d}.dat".format(snap), "rb")
      bytes = gasparts.read()
      data = struct.unpack('f'*(len(bytes)/4), bytes)
      for i in range(len(data)/6):
        maxdensity = max(maxdensity, data[6*i+3])
    print "maximal density is", maxdensity
    
    print "reading data..."
    maxtemperature = 1.
    for snap in range(snaps[0],snaps[1]+1):
      gasparts = open("gas_{0:04d}.dat".format(snap), "rb")
      bytes = gasparts.read()
      data = struct.unpack('f'*(len(bytes)/4), bytes)
      numpart = len(data)/6
      positions = [[0.,0.,0.]]*numpart
      smoothings = [0.]*numpart
      densities = [0.]*numpart
      temperatures = [0.]*numpart
      for i in range(numpart):
        positions[i] = [data[6*i], data[6*i+1], data[6*i+2]]
        smoothings[i] = data[6*i+4]
        densities[i] = data[6*i+3]
        temperatures[i] = data[6*i+5]
      
      for i in range(numpart):
        densities[i] /= maxdensity
        temperatures[i] /= maxtemperature
      glprog.add_pointset(2*snap, PointSet(positions, smoothings, densities, temperatures))
      totmemsize += glprog._pointsets[2*snap]._memsize
      
      starparts = open("stars_{0:04d}.dat".format(snap), "rb")
      bytes = starparts.read()
      data = struct.unpack('f'*(len(bytes)/4), bytes)
      numpart = len(data)/3
      positions = [[0.,0.,0.]]*numpart
      smoothings = [0.]*numpart
      densities = [0.]*numpart
      temperatures = [0.]*numpart
      for i in range(numpart):
        positions[i] = [data[3*i], data[3*i+1], data[3*i+2]]
        smoothings[i] = 0.025
        densities[i] = 0.00116
        temperatures[i] = 0.00116 # 100*log(1+temp) = 0.5 (to make stars white)
      
      glprog.add_pointset(2*snap+1, PointSet(positions, smoothings, densities, temperatures))
      totmemsize += glprog._pointsets[2*snap+1]._memsize
      
    textbox = TextBox("Temperature", gasstrength, 10.**maxTemperature, snaps[0])
    
    print "Using ", totmemsize, "bytes of GPU memory"
    inkb = totmemsize/1024
    inmb = inkb/1024
    print "This is", inkb, "KB or", inmb, "MB"
    
    glutDisplayFunc(display)
    glutReshapeFunc(resize)
    glutMouseFunc(mousepress)
    glutMotionFunc(mousemove)
    glutSpecialFunc(arrowpress)
    glutKeyboardFunc(keyfunctions)
    glutMainLoop()
    return

## Display function
#
# Wrapper around displaybuffer()
# If stereo viewing is enabled, we call displaybuffer() for the left and right buffer and
# rotate the view in between to create a 3D effect
def display():
  global view
  if view._stereo:
    glDrawBuffer(GL_BACK_LEFT)
    displaybuffer()
    view_angle = 0.034906585
#    newrmatrix = np.array([[1.,0.,0.,0.],
#                         [0.,1.,0.,0.],
#                         [0.,0.,1.,0.],
#                         [0.,0.,0.,1.]])
#    xrot = np.array([[np.cos(view_angle), 0., np.sin(view_angle), 0.],
#                     [0., 1., 0., 0.],
#                     [-np.sin(view_angle), 0., np.cos(view_angle), 0.],
#                     [0., 0., 0., 1.]])
#    newrmatrix = np.dot(newrmatrix,xrot)
#    view._rmatrix = np.dot(view._rmatrix, newrmatrix)
    # WE SHOULD WRITE A BETTER ROTATION FUNCTION FOR THIS, SINCE THIS ONE IS TOO HEAVY (WE ONLY NEED X-ROTATION)
    view.rotate_x(view_angle)
    glDrawBuffer(GL_BACK_RIGHT)
    displaybuffer()
#    deltax = -0.5
#    newrmatrix = np.array([[1.,0.,0.,0.],
#                         [0.,1.,0.,0.],
#                         [0.,0.,1.,0.],
#                         [0.,0.,0.,1.]])
#    xrot = np.array([[np.cos(-view_angle), 0., np.sin(-view_angle), 0.],
#                     [0., 1., 0., 0.],
#                     [-np.sin(-view_angle), 0., np.cos(-view_angle), 0.],
#                     [0., 0., 0., 1.]])
#    newrmatrix = np.dot(newrmatrix,xrot)
#    view._rmatrix = np.dot(view._rmatrix, newrmatrix)
    view.rotate_x(-view_angle)
  else:
    displaybuffer()
  glutSwapBuffers()

## Display function
#
# Function used to draw the contents of the window.
def displaybuffer():
    global glprog, texture, pointset, stars, startexture, gasstrength, textbox, active, maxTemperature, show_textbox, view
    # clear the contents of the view
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
    # plot gas and stars
    glUseProgram(glprog._program)
    # initialize the view
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    glTranslatef(view._center[0], view._center[1], view._center[2])
    glMultMatrixd(view._rmatrix)
    # plot gas
    glUniform1f(glprog._uniforms["strength"], 1.0)
    glUniform1f(glprog._uniforms["densityStrength"], 10.**gasstrength)
    glUniform1f(glprog._uniforms["maxTemperature"], 10.**maxTemperature)
    glActiveTexture(GL_TEXTURE0)
    glBindTexture(GL_TEXTURE_2D, texture)
    glUniform1i(glprog._uniforms["uSampler"], 0)
    glprog._pointsets[2*active].draw(glprog)
    # plot stars
    glUniform1f(glprog._uniforms["strength"], 100.)
    glUniform1f(glprog._uniforms["densityStrength"], 100.0)
    glUniform1f(glprog._uniforms["maxTemperature"], 1.)
    glActiveTexture(GL_TEXTURE0)
    glBindTexture(GL_TEXTURE_2D, startexture)
    glprog._pointsets[2*active+1].draw(glprog)
    
    if show_textbox:
      # plot the text line in the bottom left corner
      # use the default shaders
      glUseProgram(0)
      # reset the view
      glMatrixMode(GL_MODELVIEW)
      glLoadIdentity()
      # white text
      glColor3f(1., 1., 1.)
      # trial and error: position of the bottom of the box
      glRasterPos(-0.4,-0.4,-1.)
      glutBitmapString(GLUT_BITMAP_HELVETICA_12, textbox.get_text())
    return

## Function called when the window is resized
#
# We reset the viewport to the new dimensions, save the dimensions to
# the apropriate global variable and reset the perspective with the
# new aspect ratio of the window
#
# @param width New width of the window
# @param height New height of the window
def resize(width, height):
    global view
    glViewport(0, 0, width, height)
    view.resize([width, height])
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(45.0, float(width)/float(height), 0.1, 100.0)

## global variables for mouse movements
# flag to disable mousemovement actions when the left mouse button is not pressed
mousepressed = False
# position where the left mouse button was pressed, mousemovement happens relative to this
startposition = [0, 0]

## Function called when a mouse button is pressed
#
# If the left mouse button is pressed, we enable mousemovement handling and
# save the current position to the global variable
# If it is released, we disable mousemovement handling again
#
# If the mouse wheel is moved (button == 3: scroll up, 4: scroll down),
# we change the z-coordinate of the rotation center of the view
#
# @param button Integer (1-4) specifying the button that was pressed
# @param state For button 1 and 2, integer specifying if the button was pressed or released
# @param x X-coordinate of the current mouse position (in range 0-window width)
# @param y Y-coordinate of the current mouse position (in range 0-window height)
def mousepress(button, state, x, y):
    global mousepressed, startposition, view
    if button == GLUT_LEFT_BUTTON:
      if state == GLUT_DOWN:
        mousepressed = True
        startposition = [x, y]
      else:
        mousepressed = False
    if button == 3 or button == 4:
      if button == 3:
        view.translate([0., 0., 0.1])
      else:
        view.translate([0., 0., -0.1])
      glutPostRedisplay()
        
## Function called when the mouse is moved
#
# Only enabled if the left mouse button is down, as flagged by mousepress()
# We rotate the view around the x- and y-axis with an amount specified by
# the relative position of the mouse pointer w.r.t. the position where the
# left mouse button was first clicked (or where the last mousemove event occured)
# To enable continuous rotation as long as the mouse button stays down, we
# save the current mouse position to the global variable
#
# @param x X-coordinate of the current mouse position (in range 0-window width)
# @param y Y-coordinate of the current mouse position (in range 0-window height)
def mousemove(x, y):
    global mousepressed, startposition, view
    if mousepressed:
      deltax = np.pi*(x - startposition[0])/180.
      deltay = np.pi*(y - startposition[1])/180.
      startposition = [x,y]
      view.rotate(deltax, deltay)
      glutPostRedisplay()

## Function called when one of the special keys (arrows, PgUp, PgDown, End...) is pressed
#
# We only detect arrow keys
# The center of the view is moved up, down, left or right according to the arrow pressed
#
# @param key Integer specifying the key that was pressed
# @param x X-coordinate of the current mouse position (in range 0-window width)
# @param y Y-coordinate of the current mouse position (in range 0-window height)
def arrowpress(key, x, y):
  global view
  if key == GLUT_KEY_UP or key == GLUT_KEY_DOWN:
    if key == GLUT_KEY_UP:
      view.translate([0., -0.1, 0.])
    else:
      view.translate([0., 0.1, 0.])
    glutPostRedisplay()
  if key == GLUT_KEY_LEFT or key == GLUT_KEY_RIGHT:
    if key == GLUT_KEY_LEFT:
      view.translate([0.1, 0., 0.])
    else:
      view.translate([-0.1, 0., 0.])
    glutPostRedisplay()

## Function called when an ASCII-key is pressed on the keyboard
#
# Esc and some other special keys also trigger this function, see GLUT documentation
# Three distinct keys are connected to the apropriate functions
#
# @param key Integer specifying the key that was pressed
# @param x X-coordinate of the current mouse position (in range 0-window width)
# @param y Y-coordinate of the current mouse position (in range 0-window height)
def keyfunctions(key, x, y):
  if key == "e":
    # increase strength of opacity scaling in fragment shader
    change_strength(+0.1)
  if key == "d":
    # decrease strength of opacity scaling in fragment shader
    change_strength(-0.1)
  if key == "p":
    # save the current view to a .ppm image
    save_image()
  if key == "n":
    # jump to the next snapshot
    navigate(+1)
  if key == "b":
    # jump to the previous snapshot
    navigate(-1)
  if key == "w":
    # increase maximal gas temperature
    change_temperature(+0.1)
  if key == "c":
    # decrease maximal gas temperature
    change_temperature(-0.1)
  if key == "m":
    # make a movie with all snapshots and the current view state
    make_movie()

## Change the current snapshot that is displayed
#
# The amount specifies the change: +1 jumps to the next snapshot,
# -1 to the previous one.
# If the new snapshot number is out of bounds (<0 or >10), it is
# reset to the largest or smallest value to create a closed chain
#
# @param amount Amount with which the current index number that is displayed is changed
def navigate(amount):
  global active, textbox, snaps
  active += amount
  if active > snaps[1]:
    active = snaps[0]
  if active < snaps[0]:
    active = snaps[1]
  textbox._curSnap = active
  glutPostRedisplay()

## Change the maximal temperature used for the color value of the snapshot
#
# The amount specifies the change, which is in logscale (so the actual value used
# is 10**maxTemperature)
# The new value is also communicated to the TextBox
#
# @param amount Amount with which the current maximal temperature will be changed
def change_temperature(amount):
  global maxTemperature, textbox
  maxTemperature += amount
  textbox._maxTemperature = 10.**maxTemperature
  glutPostRedisplay()

## Change the strength of the opacity scaling in the fragment shader
#
# The strength is also updated in the TextBox
#
# @param amount Amount with which the strength is changed
def change_strength(amount):
  global gasstrength, textbox
  gasstrength += amount
  textbox._strength = gasstrength
  glutPostRedisplay()

## Save the current view to a .ppm image
#
# The image is written in binary format and has the same resolution as the current window
def save_image():
  global active, view
  print "Saving image"
  glPixelStorei(GL_PACK_ALIGNMENT, 1)
  data = glReadPixels(0, 0, view._dimensions[0], view._dimensions[1], GL_RGB, GL_UNSIGNED_BYTE, outputType = None)
  # for some reason, glReadPixels gives back an array with the wrong dimensions. We have to set them ourselves
  data.shape = (dimensions[1], dimensions[0], 3)
  img = open("image{0:04d}.ppm".format(active), "wb")
  img.write("P6\n" + str(view._dimensions[0]) + "\t" + str(view._dimensions[1]) + "\n255\n")
  for j in range(len(data)):
    line = data[len(data)-j-1]
    for pixel in line:
      img.write(struct.pack('B'*3, pixel[0], pixel[1], pixel[2]))
  img.close()
  print "Ready"

## Make a movie of all snapshots in the range using the current state of the view
#
# We loop over all snapshot indices and redraw the display buffer for every snapshot
# The image is then saved using save_image()
# Due to the internal workings of OpenGL, it is required for the plotting window to
# be on the foreground during this process, since otherwise the image that is saved
# will not be updated. We could remedy this by using framebuffers, but these are
# not supported by the system version of OpenGL
def make_movie():
  global active, textbox, snaps
  for i in range(snaps[0], snaps[1]+1):
    active = i
    textbox._curSnap = active
    displaybuffer()
    save_image()

# call the main function when this script is ran
if __name__ == '__main__': main()
