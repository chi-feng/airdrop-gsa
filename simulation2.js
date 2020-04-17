'use strict';

var Simulation = { };

Simulation.prototype.meta = {
  parameters: [
    { name: 'm', label: 'm', stochastic: true, distribution: { type: 'triangular', min: 4, mode: 6, max: 9, mu: 0, sigma: 1 }, default: 6 },
    { name: 'r', label: 'r', stochastic: true, distribution: { type: 'triangular', min: 0.09, mode: 0.1, max: 0.11, mu: 0, sigma: 1 }, default: 0.1 },
    { name: 'Cd', label: 'Cd', stochastic: true, distribution: { type: 'triangular', min: 0, mode: 0.5, max: 1.0, mu: 0, sigma: 1 }, default: 0 },
    { name: 'wx', label: 'wx', stochastic: true, distribution: { type: 'triangular', min: 0, mode: 0.5, max: 1.0, mu: 0, sigma: 1 }, default: 0 },
    { name: 'rp', label: 'rp', stochastic: false, default: 0.5 },
    { name: 'Cdp', label: 'Cdp', stochastic: false, default: 1.75 },
    { name: 'Fmax', label: 'Fmax', stochastic: false, default: 300 },
    { name: 'tfree', label: 'tfree', stochastic: true, distribution: { type: 'triangular', min: 0, mode: 0.5, max: 1.0, mu: 0, sigma: 1 }, default: 0 },
    { name: 'topen', label: 'topen', stochastic: true, distribution: { type: 'lognoral', min: 0, mode: 0.5, max: 1.0, mu: 0, sigma: 1 }, default: 0 },
    { name: 'x', label: 'x', stochastic: true, distribution: { type: 'normal', min: -346, mode: -340, max: -334, mu: -340, sigma: 3 }, default: 0 },
    { name: 'y', label: 'y', stochastic: true, distribution: { type: 'normal', min: 497, mode: 500, max: 503, mu: 500, sigma: 3 }, default: 0 },
    { name: 'vx', label: 'vx', stochastic: true, distribution: { type: 'normal', min: 0, mode: 0.5, max: 1.0, mu: 50, sigma: 0.5 }, default: 0 },
    { name: 'vy', label: 'vy', stochastic: true, distribution: { type: 'normal', min: 0, mode: 0.5, max: 1.0, mu: 0, sigma: 0.5 }, default: 0 },
    { name: 'xmin', label: 'xmin', stochastic: false, default: -50 },
    { name: 'xmax', label: 'xmax', stochastic: false, default: 50 },
    { name: 'nsamp', label: 'nsamp', stochastic: false, default: 1000 },
  ]
};

Simulation.prototype.getLabel = function (name) {
  this.metadata.forEach(function(group) {
    group.params.forEach(function(param, i) {
      if (name == param.name)
        return group.params[i].label;
    });
  });
};

Simulation.erf = function(x) {
  // from Chebyshev fitting formula for erf(z) from Numerical Recipes, 6.2
  var a1 =  0.254829592;
  var a2 = -0.284496736;
  var a3 =  1.421413741;
  var a4 = -1.453152027;
  var a5 =  1.061405429;
  var p  =  0.3275911;
  var t = 1.0 / (1.0 + p * Math.abs(x));
  var y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
  return Math.sign(x) * y;
};

Simulation.dudt = function(t, u, p) {
  var m     = p.m;
  var r     = p.r;
  var Cd    = p.Cd;
  var wx    = p.wx;
  var rp    = p.rp;
  var Cdp   = p.Cdp;
  var Fmax  = p.Fmax;
  var tfree = p.tfree;
  var topen = p.topen;
  var rho  = 1.22;      // atmospheric density
  var g    = 9.8;       // gravity
  var mu   = 1.81e-5;   // viscosity of air
  var x        = u[1];
  var y        = u[2];
  var vx       = u[3];
  var vy       = u[4];
  var detached = u[5];
  // calculate relative velocity
  var v = Math.sqrt((vx - wx) * (vx - wx) + vy * vy);
  // determine parachute opening status
  var parachute = 0;
  if (t < tfree || detached > 0) {
    parachute = 0;
  } else if (t > (tfree + topen)) {
    parachute = 1;
  } else {
    parachute = Math.min(1, (t - tfree) / topen);
  }
  // compute Reynolds number
  var Re = rho * v * r / mu;
  // compute atmospheric drag due to payload and parachute
  if (Math.abs(Re) > 0) {
    var Cd_eff = 24 / Re + 6 / (1 + Math.sqrt(Re)) + Cd;
    var Fd = -0.5 * rho * pi * r * r * v * v * Cd_eff;
    var Cdp_eff = 24 / Re + 6 / (1 + Math.sqrt(Re)) + Cdp;
    var Fdp = -0.5 * rho * pi * rp * rp * v * v * Cdp_eff * parachute;
  } else {
    var Fd = 0;
    var Fdp = 0;
  }
  // approximate effective drag as sum
  var ax = -Math.abs(Fd + Fdp) / m * (vx - wx) / v;
  var ay = -g - Math.abs(Fd + Fdp) / m * vy / v;
  return new Float64Array([ vx, vy, ax, ay, Math.abs(Fdp) > Fmax ? 1 : 0]);
};

Simulation.run = function(params) {
  var t = 0;
  var dt = 0.1;
  var dim = 5;
  var u = new Float64Array(dim);
  for (var i = 0; i < params.u0.length; i++) {
    u[i] = params.u0[i];
  }
  var trajectory = [ new Float64Array(u) ];
  var a = new Float64Array(dim);
  var b = new Float64Array(dim);
  var c = new Float64Array(dim);
  var d = new Float64Array(dim);
  var ub = new Float64Array(dim);
  var uc = new Float64Array(dim);
  var ud = new Float64Array(dim);
  while (true) {

    var f1 = this.dudt(t, u, params);
    for (var i = 0; i < dim; i++) {
      a[i] = dt * f1[i];
      ub[i] = u[i] + a[i] / 2;
    }
    var f2 = this.dudt(t + dt / 2, ub, params);
    for (var i = 0; i < dim; i++) {
      b[i] = dt * f2[i];
      uc[i] = u[i] + b[i] / 2;
    }
    var f3 = this.dudt(t + dt / 2, uc, params);
    for (var i = 0; i < dim; i++) {
      c[i] = dt * f3[i];
      ud[i] = u[i] + c[i];
    }
    var f4 = this.dudt(t + dt, ud, params);
    for (var i = 0; i < dim; i++) {
      d[i] = dt * f4[i];
    }
    for (var i = 0; i < dim; i++) {
      u[i] = u[i] + (a[i] * 2 * b[i] + 2 * c[i] + d[i]) / 6;
    }
    if (u[1] < 0) { // lerp to y = 0
      var x = u[0];
      var y = u[1];
      var v_x = u[2];
      var v_y = u[3];
      u[0] = x - y * v_x / v_y;
      u[1] = 0;
    }
    if (params.save_trajectory) {
      trajectory.append(new Float64Array(u));
    }
    if (u[1] == 0) {
      break;
    }
    t = t + dt;
    dt = 2 / (Math.sqrt(u[2] * u[2] + u[3] * u[3]) + 0.1);
  }
  return params.save_trajectory ? trajectory : u;
}

Simulation.getNormal = function() {
  var x, y, w;
  do {
    x = Math.random() * 2 - 1;
    y = Math.random() * 2 - 1;
    w = x * x + y * y;
  } while (w >= 1.0)
  return x * Math.sqrt(-2 * Math.log(w) / w);
};

Simulation.getLogNormal = function(mu, sigma) {
  return Math.exp(Simulation.getNormal() * Math.sqrt(sigma) + mu);
};

Simulation.getTriangular = function(a, b, c) {
  var fc = (c - a) / (b - a);
  var u = Math.random();
  if (u < fc)
    return a + Math.sqrt(u * (b - a) * (c - a));
  else
    return b - Math.sqrt((1 - u) * (b - a) * (b - c));
};

Simulation.prototype.getSamples = function() {
  var samples = { };
  this.meta.parameters.forEach(function(param) {
    if (param.stochastic ===4 true) {
      samples[param.name] = new Float64Array(params.nsamp);
      if (param.distribution.type === 'normal') {
        for (var i = 0; i < params.nsamp; i++)
          samples[param.name][i] = param.distribution.mu + param.distribution.sigma * Simulation.getNormal();
      }
      if (param.distribution.type === 'uniform') {
        for (var i = 0; i < params.nsamp; i++)
          samples[param.name][i] = Math.random() * (param.distribution.max - param.distribution.min) + param.distribution.min;
      }
      if (param.distribution.type === 'triangular') {
        for (var i = 0; i < params.nsamp; i++)
          samples[param.name][i] = Simulation.getTriangular(param.distribution.min, param.distribution.mode, param.distribution.max);
      }
      if (param.distribution.type === 'lognormal') {
        for (var i = 0; i < params.nsamp; i++)
          samples[param.name][i] = Simulation.getLogNormal(param.distribution.mu, param.distribution.sigma);
      }
    }
  });
  return samples;
};

Simulation.prototype.getDefault = function() {
  var params = { };
  this.metadata.forEach(function(paramGroup) {
    paramGroup.params.forEach(function(param) {
      params[param['name']] = param['default'];
    });
  });
  return params;
};

Simulation.prototype.parse = function(params) {
  var errors = [];
  this.metadata.forEach(function(paramGroup) {
    paramGroup.params.forEach(function(param) {
      switch (param.type) {
        case 'float': params[param.name] = parseFloat(params[param.name]); break;
        case 'int': params[param.name] = parseInt(params[param.name]); break;
      }
      if (param.min != undefined && params[param.name] < param.min) {
        errors.push(param.name + ' below minimum value of ' + param.min);
      }
      if (param.max != undefined && params[param.name] > param.max) {
        errors.push(param.name + ' above maximum value of ' + param.max);
      }
    });
  });
  params.errors = errors;
  return params;
};