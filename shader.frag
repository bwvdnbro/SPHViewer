precision mediump float;

uniform float strength;
uniform sampler2D uSampler;
uniform float densityStrength;

varying vec2 vTextureCoord;
varying float vColor;
varying float vNorm;

void main(){
  vec4 opacity = densityStrength*vNorm*texture2D(uSampler, vec2(vTextureCoord.s, vTextureCoord.t));
  float color = 100.0*log(1.0+vColor);
  gl_FragColor = vec4(strength*color, strength*(1.0-color), 1.0, opacity[0]);
}
