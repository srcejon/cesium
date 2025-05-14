import defined from "./defined.js";

/**
 * An {@link InterpolationAlgorithm} for performing Whittaker Eilers smoothing and interpolation.
 *
 * @namespace WhittakerEilersApproximation
 */
const WhittakerEilersApproximation = {
  type: "WhittakerEilers",
};

/**
 * Given the desired degree, returns the number of data points required for interpolation.
 *
 * @param {number} degree The desired degree of interpolation.
 * @returns {number} The number of required data points needed for the desired degree of interpolation.
 */
WhittakerEilersApproximation.getRequiredDataPoints = function (degree) {
  return Math.max(degree + 1.0, 2);
};

/**
 * Interpolates values using Whittaker Eilers Approximation.
 *
 * @param {number} x The independent variable for which the dependent variables will be interpolated.
 * @param {number[]} xTable The array of independent variables to use to interpolate.  The values
 * in this array must be in increasing order and the same value must not occur twice in the array.
 * @param {number[]} yTable The array of dependent variables to use to interpolate.  For a set of three
 * dependent values (p,q,w) at time 1 and time 2 this should be as follows: {p1, q1, w1, p2, q2, w2}.
 * @param {number} yStride The number of dependent variable values in yTable corresponding to
 * each independent variable value in xTable.
 * @param {number[]} [result] An existing array into which to store the result.
 * @returns {number[]} The array of interpolated values, or the result parameter if one was provided.
 */
WhittakerEilersApproximation.interpolateOrderZero = function (
  x,
  xTable,
  yTable,
  yStride,
  result,
) {

 /* console.log("x: " + x);
  for (let i = 0; i < xTable.length; i++) {
    console.log("x[" + i + "]: " + xTable[i]);
  }
  for (let i = 0; i < yTable.length; i++) {
    console.log("y[" + i + "]: " + yTable[i]);
  }
  console.log("yStride: " + yStride);*/

  if (!defined(result)) {
    result = new Array(yStride);
  }

  let i;
  let j;
  let k;
  const length = xTable.length;
  const m = xTable.length + 1;

  var xi = new Array(m);
  var yi = new Array(m);
  var w = new Array(m);
  var v1a = new Array(m);
  var v2a = new Array(m);
  var da = new Array(m * 3);
  var dtd = new Array(m * 3);
  var ca = new Array(m * 3);
  var za = new Array(m);
  var zb = new Array(m);
  var b = new Array(m);
  const lambda = 1.0;

  for (i = 0; i < yStride; i++) {
    result[i] = 0;
  }

  for (var n = 0; n < yStride; n++) {
    var insertedX = -1;
    for (i = 0, j = n, k=0; i < m; i++, j += yStride, k++) {
      xi[i] = xTable[k];
      yi[i] = yTable[j];
      w[i] = 1.0;
      if ((insertedX == -1) && (xi[k] > x)) {
        xi[i + 1] = xi[i];
        xi[i] = x;
        yi[i + 1] = yi[i];
        yi[i] = 0;
        w[i + 1] = 1.0;
        w[i] = 0.0;
        //console.log("inserted at i=" + i);
        insertedX = i;
        i++;
      }
    }

//console.log("xi: " + xi);
//console.log("yi: " + yi);
//console.log("w: " + w);
//console.log("insertedX: " + insertedX);
//console.log("m: " + m);

    for (i = 0; i < m - 1; i++) {
      v1a[i] = 1.0 / (xi[i + 1] - xi[i]);
    }
    v1a[m - 1] = 0.0;
    for (i = 0; i < m - 2; i++) {
      v2a[i] = 1.0 / (xi[i + 2] - xi[i]);
    }
    v2a[m - 1] = 0.0;
    v2a[m - 2] = 0.0;
//console.log("v1a: " + v1a);
//console.log("v2a: " + v2a);
    //Wa = w;
    // D1 = V1 * diff(I)
    for (i = 0; i < m - 1; i++) {
      da[i*3] = -v1a[i];
      da[i*3+1] = v1a[i];
    }

    // D2 = V2 * diff(D1)
    for (i = 0; i < m - 2; i++) {
      da[i*3] = v2a[i] * -da[i*3];
      da[i*3+1] = v2a[i] * (da[(i+1)*3] - da[i*3+1]);
      da[i*3+2] = v2a[i] * da[(i+1)*3+1];
    }
    for (i = 1; i <= 6; i++) {
      da[m * 3 - i] = 0;
    }

    //console.log("da: " + da);

    dtd[0 * 3] = lambda * da[0 * 3] * da[0 * 3];
    dtd[0 * 3 + 1] = lambda * da[0 * 3] * da[0 * 3 + 1];
    dtd[0 * 3 + 2] = lambda * da[0 * 3] * da[0 * 3 + 2];

    dtd[1 * 3] = lambda * (da[0 * 3 + 1] * da[0 * 3 + 1] + da[1 * 3] * da[1 * 3]);
    dtd[1 * 3 + 1] = lambda * (da[0 * 3 + 1] * da[0 * 3 + 2] + da[1 * 3] * da[1 * 3 + 1]);
    dtd[1 * 3 + 2] = lambda * (da[1 * 3] * da[1 * 3 + 2]);

    for (let row = 2; row < m; row++) {
      dtd[row * 3] = lambda * (da[(row - 2) * 3 + 2] * da[(row - 2) * 3 + 2] + da[(row - 1) * 3 + 1] * da[(row - 1) * 3 + 1] + da[row * 3] * da[row * 3]);
      dtd[row * 3 + 1] = lambda * (da[(row - 1) * 3 + 1] * da[(row - 1) * 3 + 2] + da[(row) * 3] * da[(row) * 3 + 1]);
      dtd[row * 3 + 2] = lambda * (da[(row) * 3] * da[(row) * 3 + 2]);
    }

    // Add in W
    for (i = 0; i < m; i++) {
      dtd[i * 3] = w[i] + dtd[i * 3];
    }
   // console.log("dtd: " + dtd);

    // Cholesky Decomposition

    i = 1;
    j = i - 1;
    ca[j * 3] = Math.sqrt(dtd[j * 3]);

    i = 2;
    j = i - 2;
    ca[j * 3 + 1] = 1.0 / ca[j * 3] * dtd[j * 3 + 1];

    i = 2;
    j = i - 1;
    k = i - 2;
    var sum = ca[k * 3 + 1] * ca[k * 3 + 1];
    ca[j * 3] = Math.sqrt(dtd[j * 3] - sum);

    for (i = 3; i <= m; i++) {
      j = i - 3;
      ca[j * 3 + 2] = 1.0 / ca[j * 3] * dtd[j * 3 + 2];

      k = i - 3;
      sum = ca[k * 3 + 2] * ca[k * 3 + 1];
      j = i - 2;
      ca[j * 3 + 1] = 1.0 / ca[j * 3] * (dtd[j * 3 + 1] - sum);

      k = i - 3;
      sum = ca[k * 3 + 2] * ca[k * 3 + 2];
      k = i - 2;
      sum = sum + ca[k * 3 + 1] * ca[k * 3 + 1];
      j = i - 1;
      ca[j * 3] = Math.sqrt(dtd[j * 3] - sum);
    }
    ca[m * 3 - 1] = 0;
    ca[m * 3 - 2] = 0;
    ca[m * 3 - 4] = 0;

   // console.log("ca: " + ca);

    // % Forward substitution(C' \ (w .* y))
    for (i = 0; i < m; i++) {
      b[i] = w[i] * yi[i];
    }
  //  console.log("b: " + b);


    za[0] = b[0] / ca[0];

    sum = ca[0 * 3 + 1] * za[0];
    za[1] = (b[1] - sum) / ca[1 * 3];

    for (i = 3; i <= m; i++) {
      sum = ca[(i - 3) * 3 + 2] * za[i - 3] + ca[(i - 2) * 3 + 1] * za[i - 2];
      za[i - 1] = (b[i - 1] - sum) / ca[(i - 1) * 3];
    }

 //console.log("za: " + za);

    // Backward substituion  C \ (C' \ (w .* y));
    b = za;
 //console.log("b: " + b);

    i = m - 1;
    zb[i] = b[i] / ca[i * 3];

    i = m - 2;
    sum = ca[i * 3 + 1] * zb[i + 1];
    zb[i] = (b[i] - sum) / ca[i * 3];

    for (i = m - 2; i >= 1; i--) {
      sum = ca[(i - 1) * 3 + 2] * zb[i + 1] + ca[(i - 1) * 3 + 1] * zb[i];
      zb[i-1] = (b[i - 1] - sum) / ca[(i - 1) * 3];
    }
 //console.log("zb: " + zb);


    result[n] = zb[insertedX];

    if (n==2 && result[2] > 300737642) {
      console.log("xi: " + xi);
      console.log("yi: " + yi);
      console.log("w: " + w);
      console.log("insertedX: " + insertedX);
      console.log("m: " + m);
      console.log("v1a: " + v1a);
      console.log("v2a: " + v2a);
      console.log("dtd: " + dtd);
      console.log("ca: " + ca);
      console.log("za: " + za);
      console.log("zb: " + zb);
      console.log("result[" + 2 + "]: " + result[2]);
      throw new Error();
    }

  }

  for (i = 0; i < yStride; i++) {
   // console.log("result[" + i + "]: " + result[i]);
  }

  //throw new Error();
  return result;
};
export default WhittakerEilersApproximation;
