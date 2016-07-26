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
