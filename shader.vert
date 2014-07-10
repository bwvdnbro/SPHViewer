attribute vec3 aVertexPosition;
attribute vec2 aTextureCoord;
attribute float aSmoothing;
attribute float aDensity;
attribute float aTemperature;
attribute float aDepth;

uniform float maxTemperature;

varying vec2 vTextureCoord;
varying float vColor;
varying float vNorm;

void main(){
  vec4 vertexPosition = gl_ModelViewMatrix*vec4(aVertexPosition, 1.0) + 4.0*vec4(aTextureCoord, 0.0, 1.0)*aSmoothing - 2.0*aSmoothing;
  gl_Position = gl_ProjectionMatrix * vertexPosition;
  vTextureCoord = aTextureCoord;
  vColor = aTemperature/maxTemperature;
  vNorm = 10.0*aDensity;
}
