/*
 * Copyright (C) 2023 Jacob Farnsworth
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

const { quat, vec2, vec3, vec4, mat3, mat4 } = glMatrix;

/**
 * Modulus function that outputs only positive values
 * @param {*} n 
 * @param {*} m 
 * @returns n modulo m
 */
function mod(n, m) {
  return ((n % m) + m) % m;
}

/**
 * Check if two numbers are equal within a certain
 * tolerance range. By default, 0.01.
 * @param {*} left the left number to check
 * @param {*} right the right number to check
 * @returns whether the numbers are equal
 */
function tolerantEquals(left, right) {
  return Math.abs(left - right) < 0.01;
}

/**
 * Extracts everything after the first character in a string
 * as a number.
 * @param {*} str String to extract digits from
 * @returns the number
 */
function takeDigits(str) {
  // Check if the string is empty or has less than 2 characters
  if (!str || str.length < 2) {
    return NaN;
  }

  // Extract the substring from the second character to the end
  var numberPart = str.slice(1);

  // Convert the extracted substring to a number
  return Number(numberPart);
}

/**
 * Compute the coordinates of all the
 * complex nth roots of unity.
 * @param {Number} n 
 * @returns array containing pairs of coordinates
 */
function nthRootsOfUnity(n) {
  let roots = [];
  for (let k = 0; k < n; k++) {
      let angle = 2 * Math.PI * k / n;
      let x = Math.cos(angle);
      let y = Math.sin(angle);
      roots.push([x, y]);
  }
  return roots;
}

/**
 * Clamps a value between `min` and `max`.
 * @param {Number} number the value to clamp
 * @param {Number} min minimum value
 * @param {Number} max maximum value
 * @returns the clamped value
 */
function clamp(number, min, max) {
  return Math.max(min, Math.min(number, max));
}

/**
 * Convert degrees to radians.
 * @param {Number} deg the angle in degrees
 * @returns the angle in radians
 */
function deg2rad(deg) {
  return Math.PI * deg / 180.0;
}

/**
 * Convert radians to degrees.
 * @param {Number} rad the angle in radians
 * @returns the angle in degrees
 */
function rad2deg(rad) {
  return 180.0 / Math.PI * rad;
}

/**
 * Creates a vertex color array for a ball with specified
 * radius and color.
 * @param {*} radius Radius of the ball
 * @param {*} color Color of the ball
 * @returns the vertex color array
 */
function createBallVertexArray(radius, color) {
  var vertColorArray = [0.0, 0.0, 0.0, color[0], color[1], color[2]];

  const numVertsInBall = 16;
  let angle = 0.0;
  const angleIncrement = Math.PI * 2.0 / numVertsInBall;

  // Add a vertex along the edge for each
  // angle increment
  for(let k = 0; k <= numVertsInBall; ++k) {
    const ballVertX = radius * Math.cos(angle);
    const ballVertY = radius * Math.sin(angle);
    vertColorArray.push(ballVertX, ballVertY, 0.0, color[0], color[1], color[2]);

    angle += angleIncrement;
  }

  return vertColorArray;
}

// Vertices and Colors (R, G, B)
var polygonVerticesColors = new Float32Array([]);

const ballVerticesColors = [
  new Float32Array(createBallVertexArray(0.06, [1.0, 0.0, 0.0])),
  new Float32Array(createBallVertexArray(0.06, [0.0, 1.0, 0.0])),
  new Float32Array(createBallVertexArray(0.06, [0.0, 0.0, 1.0])),
  new Float32Array(createBallVertexArray(0.06, [1.0, 1.0, 0.0])),
  new Float32Array(createBallVertexArray(0.06, [0.0, 1.0, 1.0])),
  new Float32Array(createBallVertexArray(0.06, [1.0, 0.0, 1.0])),
  new Float32Array(createBallVertexArray(0.06, [0.78, 0.87, 0.98])),
  new Float32Array(createBallVertexArray(0.06, [0.57, 0.93, 0.07])),
  new Float32Array(createBallVertexArray(0.06, [0.94, 0.52, 0.41])),
  new Float32Array(createBallVertexArray(0.06, [0.57, 0.02, 0.62])),
  new Float32Array(createBallVertexArray(0.06, [0.80, 0.46, 0.78])),
  new Float32Array(createBallVertexArray(0.06, [0.68, 0.36, 0.44])),
];

var numPolygonVertices = 5; // Number of vertices
var numPolygonSides = 3;

var currentActionProgress = 0.0;
var currentActionSpeed = 0.0;
var currentActionMatrix = mat4.create();
var currentStateMatrix = mat4.create();

var targetAngle = 0.0;
var targetReflectionAngle = 0.0;

/**
 * Compute a rotation transformation matrix rotating counterclockwise
 * `targetAngle` degrees.
 * @param {*} targetAngle the rotation angle in degrees
 * @returns the transformation matrix
 */
function computeTargetRotationMatrix(targetAngle) {
  const cosTh = Math.cos(deg2rad(targetAngle));
  const sinTh = Math.sin(deg2rad(targetAngle));

  const rotationMatrix = mat4.fromValues(
    cosTh, sinTh, 0.0, 0.0,
    -sinTh, cosTh, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0
  );

  return rotationMatrix;
}

/**
 * Compute a reflection transformation matrix across
 * the line extending from the origin at `reflectionAngle` degrees.
 * @param {*} reflectionAngle the angle of the reflection axis
 * @returns the transformation matrix
 */
function computeTargetReflectionMatrix(reflectionAngle) {
  const cosTh = Math.cos(deg2rad(reflectionAngle * 2.0));
  const sinTh = Math.sin(deg2rad(reflectionAngle * 2.0));

  const reflectionMatrix = mat4.fromValues(
    cosTh, sinTh, 0.0, 0.0,
    sinTh, -cosTh, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0
  );

  return reflectionMatrix;
}

/**
 * Interpolates a rotation matrix between 0 and `rotationAngles` degrees
 * according to `progress`, which ranges between 0 and 1. 
 * @param {*} rotationAngle the rotation angle in degrees
 * @param {*} progress progress of the action between 0 and 1
 * @returns the rotation matrix
 */
function interpolateRotationMatrix(rotationAngle, progress) {
  var rotationMatrix = computeTargetRotationMatrix(rotationAngle * progress);
  return rotationMatrix;
}

/**
 * Interpolates a reflection matrix across the line extending from the origin
 * at `axisAngle` degrees, according to `progress1`, which ranges between 0 and 1.
 * @param {*} reflectionAxisAngle the reflection axis angle in degrees
 * @param {*} progress progress of the action between 0 and 1
 * @returns the reflection transformation matrix
 */
function interpolateReflectionMatrix(reflectionAxisAngle, progress) {
  var reflectionMatrix = computeTargetReflectionMatrix(reflectionAxisAngle);
  var resultMatrix = mat4.create();

  var scaledReflectionMatrix = mat4.clone(reflectionMatrix);
  var scaledIdentityMatrix = mat4.create();

  mat4.multiplyScalar(scaledReflectionMatrix, reflectionMatrix, progress);
  mat4.multiplyScalar(scaledIdentityMatrix, scaledIdentityMatrix, 1.0 - progress);

  mat4.add(resultMatrix, scaledReflectionMatrix, scaledIdentityMatrix);

  return resultMatrix;
}

/**
 * Interpret the current state as either a rotation
 * or a reflection and get the angle of rotation and reflection
 * respectively.
 * @returns data about the current state
 */
function interpretCurrentState() {
  const m00 = currentStateMatrix[0];
  const m01 = currentStateMatrix[1];
  const m10 = currentStateMatrix[4];
  const m11 = currentStateMatrix[5];

  var identityMatrix = mat4.create();
  var stateTransposeMatrix = mat4.create();
  var stateInverseMatrix = mat4.create();
  mat4.transpose(stateTransposeMatrix, currentStateMatrix);
  mat4.invert(stateInverseMatrix, currentStateMatrix);

  let reflection = false;
  let rotation = false;
  let reflectionAngle = 0;
  let rotationAngle = 0;

  // Compute eigenvalues of the matrix
  // First get coefficients of the quadratic
  // to solve
  const a = 1.0;
  const b = -m00 - m11;
  const c = m00 * m11 - m01 * m10;

  // Compute eigenvalues
  const sqrtdet = Math.sqrt(b * b - 4.0 * a * c)
  const x0 = (-b + sqrtdet) / (2.0 * a);
  const x1 = (-b - sqrtdet) / (2.0 * a);

  let eigenvectors = [];

  // Helper function to compute an eigenvector for a given eigenvalue
  function computeEigenvectorForEigenvalue(lambda) {
    let a = m00 - lambda;
    let b = m01;
    let c = m10;
    let d = m11 - lambda;

    // Check for the (1, 0) special case
    if (a === 0 && c === 0) {
        return [1, 0];
    }

    // Otherwise, solve for the second component assuming the first component is 1
    let v2 = NaN;
    if (b !== 0) {
        v2 = -a / b;
    } else if (d !== 0) {
        v2 = -c / d;
    } else {
        // do nothing
    }
    return [1, v2];
  }

  // Compute eigenvectors for each eigenvalue
  eigenvectors.push(computeEigenvectorForEigenvalue(x0));
  eigenvectors.push(computeEigenvectorForEigenvalue(x1));
  let eigenvector;
  let eigenvalue;

  // If the matrix represents a reflection,
  // one of the eigenvalues will be -1, and its corresponding
  // eigenvector will have the axis of rotation
  if(tolerantEquals(x0, x1)) {
    // do nothing, proceed to check rotation later
    // (can end up here if the transformation is a rotation
    // of even order)
  } else if(tolerantEquals(x0, -1.0)) {
    reflection = true;
    eigenvector = eigenvectors[1];
  } else if(tolerantEquals(x1, -1.0)) {
    reflection = true;
    eigenvector = eigenvectors[0];
  } else {
    // do nothing
  }

  if(reflection) {
    // Compute the angle from the origin
    reflectionAngle = rad2deg(Math.atan2(eigenvector[1], eigenvector[0]));
  } else {
    // In this case, we have a rotation, in which case we can
    // get the angle by applying the rotation to a temporary vector,
    // taking the dot product with the original vector to obtain the angle.

    rotation = true;
    const unitVec = vec4.fromValues(1.0, 0.0, 0.0, 0.0);
    var tmp = vec4.create();
    vec4.transformMat4(tmp, unitVec, currentStateMatrix);

    // vectors (should be) already normalized
    rotationAngle = rad2deg(Math.atan2(tmp[1], tmp[0]));
  }

  return { reflection, rotation, reflectionAngle, rotationAngle };
}

/**
 * Get a string representing the current group element by
 * interpreting the current state matrix.
 * @returns string corresponding to the current element
 */
function getElementString() {
  const { reflection, rotation, reflectionAngle, rotationAngle } = interpretCurrentState();

  let elem;

  if(reflection) {
    elem = "s";
    const digit = mod(Math.round(reflectionAngle * 2.0 / (360.0 / numPolygonSides)), numPolygonSides);
    elem += digit;
  } else {
    elem = "r";
    const digit = mod(Math.round(rotationAngle / (360.0 / numPolygonSides)), numPolygonSides);
    elem += digit;
  }

  return elem;
}

function isAnimationInProgress() {
  return currentActionSpeed > 0.0;
}

function initMatrices(gl) {
  // Perspective Projection Matrix
  var projectionMatrix = mat4.create();
  mat4.perspective(projectionMatrix, 45 * Math.PI / 180, gl.canvas.width / gl.canvas.height, 0.1, 100);

  // View Matrix (Camera)
  var viewMatrix = mat4.create();
  mat4.lookAt(viewMatrix, [0, 0, 3], [0, 0, 0], [0, 1, 0]);

  return { projectionMatrix, viewMatrix };
}

function init() {
  // Initialize the triangle group
  groupNew();

  // look up the elements we want to affect
  var stateElement = document.querySelector("#elemState");
  var stateDescElement = document.querySelector("#elemDesc");
  
  // Create text nodes to save some time for the browser.
  var stateNode = document.createTextNode("");
  var stateDescNode = document.createTextNode("");
  
  // Add those text nodes where they need to go
  stateElement.appendChild(stateNode);
  stateDescElement.appendChild(stateDescNode);

  // Get canvas object from the DOM
  var canvas = document.getElementById("glcanvas");

  // Init WebGL context
  var gl = canvas.getContext("webgl");
  if (!gl) {
      console.log("Failed to get the rendering context for WebGL");
      return;
  }

  // Init shaders
  var vs = document.getElementById('shaderVs').innerHTML;
  var fs = document.getElementById('shaderFs').innerHTML;
  if (!initShaders(vs, fs)) {
      console.log('Failed to intialize shaders.');
      return;
  }

  // Create a buffer object
  var vertexColorBuffer = gl.createBuffer();
  if (!vertexColorBuffer) {
      console.log('Failed to create the buffer object');
      return -1;
  }

  var matrices = initMatrices(gl);

  drawScene();

  function drawScene() {
    stateDescNode.nodeValue = getElementString();

    var modelMatrix = mat4.create();

    setMatrixUniforms(matrices.projectionMatrix, matrices.viewMatrix);

    updateAnimationValues();

    // Clear canvas
    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT);
  
    // Render the polygon
    setModelMatrix(modelMatrix);
    renderTriangles(polygonVerticesColors, numPolygonVertices);

    // Render vertex points of the polygon
    renderPolygonPoints();

    // Call drawScene again next frame
    requestAnimationFrame(drawScene);
  }

  function renderPolygonPoints() {
    for(let k = 0; k < numPolygonSides; ++k) {
      const theta = k * 2.0 * Math.PI / numPolygonSides;
      const ptX = Math.cos(theta);
      const ptY = Math.sin(theta);

      const translationVector = vec4.fromValues(ptX, ptY, 0.0, 0.0);
      var modelMatrix = mat4.create();
      mat4.fromTranslation(modelMatrix, translationVector);

      setModelMatrix(modelMatrix);
      renderBall(k);
    }
  }

  function resetAnimations() {
    currentActionProgress = 0.0;
    currentActionSpeed = 0.0;
    currentActionMatrix = mat4.create();
    targetAngle = 0.0;
    targetReflectionAngle = 0.0;
  }

  function updateAnimationValues() {
    if(currentActionSpeed > 0.0) {
      currentActionProgress = clamp(currentActionProgress + currentActionSpeed, 0.0, 1.0);

      if(targetAngle > 0.0) {
        // Progress a rotation action

        currentActionMatrix = interpolateRotationMatrix(targetAngle, currentActionProgress);
      } else {
        currentActionMatrix = interpolateReflectionMatrix(targetReflectionAngle, currentActionProgress);
      }

      if(currentActionProgress == 1.0) {
        // We have reached the end of the action, update the state matrix
        mat4.multiply(currentStateMatrix, currentActionMatrix, currentStateMatrix);

        resetAnimations();
      }
    }

    // If animation is in progress, we want to
    // hide the table buttons
    if(isAnimationInProgress()) {
      document.getElementById("tblGrpActButtons").style.visibility = "hidden";
    } else {
      document.getElementById("tblGrpActButtons").style.visibility = "";
    }
  }

  /**
   * Set the projection and view matrix uniforms.
   * @param {*} projectionMatrix the projection matrix
   * @param {*} viewMatrix the view matrix
   */
  function setMatrixUniforms(projectionMatrix, viewMatrix) {
    var u_Projection = gl.getUniformLocation(gl.program, 'u_Projection');
    var u_View = gl.getUniformLocation(gl.program, 'u_View');
    gl.uniformMatrix4fv(u_Projection, false, projectionMatrix);
    gl.uniformMatrix4fv(u_View, false, viewMatrix);
    
  }

  /**
   * Set the model matrix uniform.
   * @param {*} modelMatrix the model matrix
   */
  function setModelMatrix(modelMatrix) {
    var u_Model = gl.getUniformLocation(gl.program, 'u_Model');
    gl.uniformMatrix4fv(u_Model, false, modelMatrix);
  }

  function renderBall(ballIndex) {
    var verticesColors = ballVerticesColors[ballIndex % ballVerticesColors.length];

    renderTriangles(verticesColors, verticesColors.length / 6);
  }

  /**
   * Updates the vertex buffers and renders
   * with TRIANGLE_FAN
   * @param {*} verticesColors Vertex color array
   * @param {*} numVertices number of vertex-color pairs in the array
   * @returns 
   */
  function renderTriangles(verticesColors, numVertices) {
    // Bind the buffer object to target
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexColorBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, verticesColors, gl.STATIC_DRAW);

    // Write the positions of vertices to a vertex shader
    var n = updateVertexBuffers(verticesColors);
    if (n < 0) {
        console.log('Failed to set the positions of the vertices');
        return;
    }

    // Draw
    gl.drawArrays(gl.TRIANGLE_FAN, 0, numVertices);
  }
  
  /**
   * Update shader uniforms and set the render transform matrix.
   * @param {*} verticesColors vertex color array
   * @returns negative if failed, positive otherwise
   */
  function updateVertexBuffers(verticesColors) {
    var FSIZE = verticesColors.BYTES_PER_ELEMENT;

    // Assign the buffer object to a_Position and enable the assignment
    var a_Position = gl.getAttribLocation(gl.program, 'a_Position');
    if(a_Position < 0) {
        console.log('Failed to get the storage location of a_Position');
        return -1;
    }
    gl.vertexAttribPointer(a_Position, 3, gl.FLOAT, false, FSIZE * 6, 0);
    gl.enableVertexAttribArray(a_Position);

    // Assign the buffer object to a_Color and enable the assignment
    var a_Color = gl.getAttribLocation(gl.program, 'a_Color');
    if (a_Color < 0) {
        console.log('Failed to get the storage location of a_Color');
        return -1;
    }
    gl.vertexAttribPointer(a_Color, 3, gl.FLOAT, false, FSIZE * 6, FSIZE * 3);
    gl.enableVertexAttribArray(a_Color);

    // Obtain the actual transformation matrix
    // by multiplying the current action matrix by the state matrix
    var actualTransformationMatrix = mat4.create();
    mat4.multiply(actualTransformationMatrix, currentActionMatrix, currentStateMatrix);

    // Set rotation uniform
    var u_Transform = gl.getUniformLocation(gl.program, 'u_Transform');
    gl.uniformMatrix4fv(u_Transform, false, actualTransformationMatrix);

    return 1;
  }
  
  function initShaders(vs_source, fs_source) {
    // Compile shaders
    var vertexShader = makeShader(vs_source, gl.VERTEX_SHADER);
    var fragmentShader = makeShader(fs_source, gl.FRAGMENT_SHADER);
  
    // Create program
    var glProgram = gl.createProgram();
  
    // Attach and link shaders to the program
    gl.attachShader(glProgram, vertexShader);
    gl.attachShader(glProgram, fragmentShader);
    gl.linkProgram(glProgram);
    if (!gl.getProgramParameter(glProgram, gl.LINK_STATUS)) {
        alert("Unable to initialize the shader program");
        return false;
    }
  
    // Use program
    gl.useProgram(glProgram);
    gl.program = glProgram;
  
    return true;
  }
  
  function makeShader(src, type) {
    var shader = gl.createShader(type);
    gl.shaderSource(shader, src);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        alert("Error compiling shader: " + gl.getShaderInfoLog(shader));
        return;
    }
    return shader;
  }
}

/**
 * Initializes vertex and color array data
 * according to number of sides.
 * @param {*} numSides the number of sides of the polygon
 */
function initializePolygon(numSides) {
  var nthRootCoords = nthRootsOfUnity(numSides);
  var colors = [
    [0.7, 0.7, 0.7],
    [0.3, 0.3, 0.3]
  ];

  var vertColorArray = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
  var i = 0;
  for(const root of nthRootCoords) {
    const color = colors[(i++) % colors.length];
    vertColorArray.push(root[0], root[1], 0, color[0], color[1], color[2]);
  }

  // Add a copy of the first vertex
  vertColorArray.push(
    vertColorArray[6], vertColorArray[7], vertColorArray[8],
    vertColorArray[9], vertColorArray[10], vertColorArray[11]);

  // change global verticesColors to the new array
  polygonVerticesColors = new Float32Array(vertColorArray);

  // change global numVertices to the new
  numPolygonSides = numSides;
  numPolygonVertices = nthRootCoords.length + 2;
}

function populateActionTable(numSides) {
  // Populate the rotations column with buttons
  var rotationsColumn = document.getElementById("tblGrpColRotations");
  rotationsColumn.innerHTML = "";
  for(let k = 0; k < numSides; ++k) {
    rotationsColumn.innerHTML += `<button onclick="groupAction('r${k}')">r${k}</button>`;
  }

  // Populate the reflections column with buttons
  var reflectionsColumn = document.getElementById("tblGrpColReflections");
  reflectionsColumn.innerHTML = "";
  for(let k = 0; k < numSides; ++k) {
    reflectionsColumn.innerHTML += `<button onclick="groupAction('s${k}')">s${k}</button>`;
  }
}

function groupStateReset() {
  console.log(`Resetting state`);

  // Reset the state matrix
  currentStateMatrix = mat4.create();
}

/**
 * Triggered when the user chooses a group action
 * @param {*} action 
 */
function groupAction(action) {
  console.log(`Applying action: (${action})`);

  if(action === 'r0') {
    // identity action: do nothing
    return;
  }

  const actionIndex = takeDigits(action);

  if(action[0] === 's') {
    currentActionSpeed = 0.01667;
    targetReflectionAngle = 1.0 / 2.0 * actionIndex * 360.0 / numPolygonSides;
  } else if(action[0] === 'r') {
    currentActionSpeed = 0.01667;
    targetAngle = actionIndex * 360.0 / numPolygonSides;
  } else {
    console.log(`Unrecognized group action ${action}`);
  }
}

/**
 * Triggered when the user selects # of sides in the dropdown
 */
function groupNew() {
  var selector = document.getElementById("polySelector");
  const numSides = Number(selector.value);

  // Initialize the vertex buffer, etc based on the number
  // of sides
  initializePolygon(numSides);

  // Reset the group state so that we start with the identity
  // state
  groupStateReset();

  // Reset the group actions table and fill it
  // with buttons for all the rotation and reflection
  // actions
  populateActionTable(numSides);
}