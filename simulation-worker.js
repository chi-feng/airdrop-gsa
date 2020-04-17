importScripts('simulation.js');

function run_mc(event) {
  var sim = new Simulation();
  var params = sim.parse(event.data);
  if (params.errors.length > 0) {
    postMessage({action:'error', errors: params.errors});
    return;
  }
  var x = [], v = [], p = [], b = [];
  var trajectories = [];
  var samples = sim.getSamples(params);
  for (var i = 0; i < params.n_samp; i++) {

    params.w_x     = samples.w_x[i];
    params.t_open  = samples.t_open[i];
    params.C_d     = samples.C_d[i];
    params.m       = samples.m[i];
    params.v_x     = samples.v_x[i];

    params.save_trajectory = (i < 25);
    var result = sim.run(params);
    if (params.save_trajectory) {
      if (i == 1)
        console.log(result)
      var traj_x = [];
      var traj_y = [];
      var traj_open_x = [];
      var traj_open_y = [];
      var open_x, open_y = -1;
      var full_x, full_y = -1;
      var land_x, land_y = -1;
      var torn_x, torn_y = -1;
      var success = 0;
      for (var j = 0; j < result.length; j++) {
        if (result[j].p == 0) {
          traj_x.push(result[j].x);
          traj_y.push(result[j].y);
        } else {
          if (traj_open_x.length == 0) {
            traj_x.push(result[j].x);
            traj_y.push(result[j].y);
          }
          traj_open_x.push(result[j].x);
          traj_open_y.push(result[j].y);
        }
        if (result[j].p > 0 && open_y == -1) {
          open_x = result[j].x;
          open_y = result[j].y;
        }
        if (result[j].p == 1 && full_y == -1) {
          full_x = result[j].x;
          full_y = result[j].y;
        }
        if (result[j].b > 0 && torn_y == -1) {
          torn_x = result[j].x;
          torn_y = result[j].y;
        }
      }
      var last = result[result.length - 1];
      var land_x = last.x, land_y = last.y;
      var land_v = Math.sqrt(last.v_x * last.v_x + last.v_y * last.v_y);
      if (land_v <= params.v_max)
        success = 1;
      trajectories.push({
        before:{x: traj_x, y: traj_y},
        after:{x:traj_open_x, y:traj_open_y},
        open:{x:open_x,y:open_y},
        full:{x:full_x,y:full_y},
        torn:{x:torn_x,y:torn_y},
        land:{x:land_x,y:land_y,v:land_v},
        success: success
      });
    }
    var last = result[result.length - 1];
    x.push(last.x);
    v.push(Math.sqrt(last.v_x * last.v_x + last.v_y * last.v_y));
    p.push(last.p);
    b.push(last.b);
  }
  postMessage({action:'run_mc', x: x, v: v, p: p, b:b, params: params, samples: samples, trajectories: trajectories});
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
    var y_A = {x: new Float64Array(N), v: new Float64Array(N)};
    var y_B = {x: new Float64Array(N), v: new Float64Array(N)};
    var y_C = {x: new Float64Array(N), v: new Float64Array(N)};
    for (var i = 0; i < N; i++) {
      for (var j = 0; j < keys.length; j++) {
        params[keys[j]] = A[i * keys.length + j];
      }
      var trajectory = sim.run(params);
      y_A.x[i] = trajectory[trajectory.length-1].x;
      y_A.v[i] = Math.sqrt(Math.pow(trajectory[trajectory.length-1].v_x,2)+Math.pow(trajectory[trajectory.length-1].v_y,2));
      for (var j = 0; j < keys.length; j++) {
        params[keys[j]] = B[i * keys.length + j];
      }
      var trajectory = sim.run(params);
      y_B.x[i] = trajectory[trajectory.length-1].x;
      y_B.v[i] = Math.sqrt(Math.pow(trajectory[trajectory.length-1].v_x,2)+Math.pow(trajectory[trajectory.length-1].v_y,2));
    }
    var f0 = {x:0, v:0};
    for (var i = 0; i < N; i++) {
      f0.x += y_A.x[i];
      f0.v += y_A.v[i];
    }
    f0.x = f0.x / N;
    f0.v = f0.v / N;
    // main and total effect sensitivities
    var S = {x: new Float64Array(keys.length), v: new Float64Array(keys.length)};
    var ST = {x: new Float64Array(keys.length), v: new Float64Array(keys.length)};

    function mean(arr) {
      var sum = 0;
      for (var i = 0; i < arr.length; i++) {
        sum += arr[i];
      }
      return sum / arr.length;
    };

    function variance(arr) {
      var sum = 0;
      var mu = mean(arr);
      for (var i = 0; i < arr.length; i++) {
        sum += (arr[i] - mu) * (arr[i] - mu);
      }
      return sum / (arr.length - 1);
    };

    var muA = {x:mean(y_A.x), v:mean(y_A.v)};
    var muB = {x:mean(y_B.x), v:mean(y_B.v)};
    var varA = {x:variance(y_A.x), v:variance(y_A.v)};
    var varB = {x:variance(y_B.x), v:variance(y_B.v)};
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
        y_C.x[i] = trajectory[trajectory.length-1].x;
        y_C.v[i] = Math.sqrt(Math.pow(trajectory[trajectory.length-1].v_x,2)+Math.pow(trajectory[trajectory.length-1].v_y,2));
      }

      console.log(k, keys[k]);

      S.x[k] = (2 * N / (2 * N - 1) * (dot(y_A.x, y_C.x) / N - Math.pow(muA.x + muB.x, 2) / 4 + (varA.x + varB.x) / (4 * N))) / (varA.x);
      var sumsq = 0; for (var i = 0; i < N; i++) { sumsq += Math.pow(y_B.x[i] - y_C.x[i], 2); }
      ST.x[k] = 1 / (2*N) * sumsq / varA.x;

      S.v[k] = (2 * N / (2 * N - 1) * (dot(y_A.v, y_C.v) / N - Math.pow(muA.v + muB.v, 2) / 4 + (varA.v + varB.v) / (4 * N))) / (varA.v);
      var sumsq = 0; for (var i = 0; i < N; i++) { sumsq += Math.pow(y_B.v[i] - y_C.v[i], 2); }
      ST.v[k] = 1 / (2*N) * sumsq / varA.v;

      /*
      S.x[k] = (dot(y_A.x, y_C.x) / N - f0.x * f0.x) / (dot(y_A.x, y_A.x) / N - f0.x * f0.x);
      ST.x[k] = 1 - (dot(y_B.x, y_C.x) / N - f0.x * f0.x) / (dot(y_A.x, y_A.x) / N - f0.x * f0.x);

      S.v[k] = (dot(y_A.v, y_C.v) / N - f0.v * f0.v) / (dot(y_A.v, y_A.v) / N - f0.v * f0.v);
      ST.v[k] = 1 - (dot(y_B.v, y_C.v) / N - f0.v * f0.v) / (dot(y_A.v, y_A.v) / N - f0.v * f0.v);
      */
    }
    return { 'keys': keys, 'main_effect_sensitivites': S, 'total_effect_sensitivites': ST };
  }

  var sim = new Simulation();
  var params = sim.parse(event.data);

  if (params.errors.length > 0) {
    postMessage({action:'error', errors: params.errors});
    return;
  }

  var result = compute_sensitivities(params)
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