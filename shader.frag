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
