importScripts('simulation2.js');

function run_mc(event) {
  var sim = new Simulation();

  var params = sim.getParams(event.data);
  if (params.errors.length > 0) {
    postMessage({action:'error', errors: params.errors});
    return;
  }

  var x = [], v = [], p = [], b = [];
  var trajectories = [];
  var samples = sim.getSamples(params);

  for (var i = 0; i < params.nsamp; i++) {

    for (var j = 0; j < sim.meta.parameters.length; j++) {
      var param = sim.meta.parameters[j];
      if (param.stochastic === true) {
        params[param.name] = samples[param.name][i];
      }
    }

    if (i < 25) params.save_trajectory = true;
    var result = sim.run(params);

    if (params.save_trajectory) {
      var x = [], y = [];
      for (var j = 0; j < result.length; j++) {
        x.push(result[j][0]);
        y.push(result[j][1]);
      }
      trajectories.push({
        complete: {x: x, y: y}
      });
      result = result[result.length - 1];
    }
    xf.push(result[0]);
    vf.push(Math.sqrt(result[2] * result[2] + result[3] * result[3]));

    var mu = 2.3;
    var sigma = 0.15;
    var thresh = 0.5 + 0.5 * Simulation.erf((Math.log(vf[vf.length - 1]) - mu) / (Math.sqrt(2) * sigma));
    if (Math.random() < thresh)
      intact.push(0);
    else
      intact.push(1);
  }
  postMessage({action:'run_mc', xf: xf, vf: vf, intact: intact, params: params, samples: samples, trajectories: trajectories});
}

function run_sa(event) {

  var dot = function (a, b) {
    var prod = 0;
    for (var i = 0; i < a.length; i++)
      prod += a[i] * b[i];
    return prod;
  };

  var compute_sensitivities = function(params) {
    var N = params.n_samp;
    A_samps = sim.getSamples(params);
    var keys = Object.keys(A_samps);
    B_samps = sim.getSamples(params);
    var A = new Float64Array(keys.length * N);
    var B = new Float64Array(keys.length * N);
    keys.forEach(function(key, j) {
      var A_samp = A_samps[key];
      var B_samp = B_samps[key];
      for (var i = 0; i < N; i++) {
        A[i * keys.length + j] = A_samp[i];
        B[i * keys.length + j] = B_samp[i];
      }
    });
    var y_A = new Float64Array(N);
    var y_B = new Float64Array(N);
    var y_C = new Float64Array(N);
    for (var i = 0; i < N; i++) {
      for (var j = 0; j < keys.length; j++) {
        params[keys[j]] = A[i * keys.length + j];
      }
      var trajectory = sim.run(params);
      y_A[i] = trajectory[trajectory.length-1].x;
      for (var j = 0; j < keys.length; j++) {
        params[keys[j]] = B[i * keys.length + j];
      }
      var trajectory = sim.run(params);
      y_B[i] = trajectory[trajectory.length-1].x;
    }
    var f0 = 0;
    for (var i = 0; i < N; i++)
      f0 += y_A[i];
    f0 = f0 / N;
    // main and total effect sensitivities
    var S = new Float64Array(keys.length);
    var ST = new Float64Array(keys.length);
    for (var k = 0; k < keys.length; k++) {
      var C = Float64Array.from(B);
      for (var i = 0; i < N; i++) {
        C[i * keys.length + k] = A[i * keys.length + k];
      }
      for (var i = 0; i < N; i++) {
        for (var j = 0; j < keys.length; j++) {
          params[keys[j]] = C[i * keys.length + j];
        }
        var trajectory = sim.run(params);
        y_C[i] = trajectory[trajectory.length-1].x;
      }
      S[k] = (dot(y_A, y_C) / N - f0 * f0) / (dot(y_A, y_A) / N - f0 * f0);
      ST[k] = 1 - (dot(y_B, y_C) / N - f0 * f0) / (dot(y_A, y_A) / N - f0 * f0);
    }
    return { 'keys': keys, 'main_effect_sensitivites': S, 'total_effect_sensitivites': ST };
  }

  var sim = new Simulation();
  var params = sim.parse(event.data);

  if (params.errors.length > 0) {
    postMessage({action:'error', errors: params.errors});
    return;
  }

  var result = compute_sensitivities(params);
  result.action = 'run_sa';
  postMessage(result);

}

onmessage = function(event) {

  if (event.data.action == 'run_mc') {
    var t0 = performance.now();
    run_mc(event);
    console.log("run_mc() took " + (performance.now() - t0) + " ms.");
  }

  if (event.data.action == 'run_sa') {
    var t0 = performance.now();
    run_sa(event);
    console.log("run_sa() took " + (performance.now() - t0) + " ms.");
  }

};