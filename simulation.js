'use strict';

function Simulation() { }

Simulation.prototype.metadata = [
  {
    label: 'Payload Design', params: [
      { name: 'r', type: 'float', label: 'r', default: 0.1, units: 'm', readonly: true},
      { name: 'C_d', type: 'float', label: 'C<sub>D</sub>', default: 0.5, nobreak: true, min: 0.2, max: 0.7, readonly: true},
      { name: 'C_derr', type: 'float', label: '&plusmn;', default: 0.1, min: 0.05, max: 0.2 },
      { name: 'mean_m', type: 'float', label: 'm', default: 5, nobreak: true, min: 4, max: 9, readonly: true},
      { name: 'var_m', type: 'float', label: '&plusmn;', default: 0.5, max: 1.0, min: 0.25},
    ],
  }, {
    label: 'Parachute Design', params: [
      { name: 'r_p', type: 'float', label: 'r<sub>p</sub>', default: 0.5, units: 'm', min: 0.25, max: 1.0, readonly: true},
      { name: 'C_dp', type: 'float', label: 'C<sub>Dp</sub>', default: 1.75, readonly: true},
      { name: 't_p', type: 'float', label: 't<sub>freefall</sub>', default: 9, units: 's', min: 0, max: 20},
      { name: 'mean_t_open', type: 'float', label: 't<sub>open</sub>', default: 4.5, nobreak: true, min: 1, max: 10},
      { name: 'var_t_open', type: 'float', label: '&plusmn;', default: 2, min: 1, max: 4.0},
    ]
  }, {
    label: 'Initial Conditions', params: [
      { name: 'x', type: 'float', label: 'x', default: -320, units: 'm' },
      { name: 'y', type: 'float', label: 'y', default: 500, units: 'm' , readonly: true},
      { name: 'mean_v_x', type: 'float', label: 'v<sub>x</sub>', default: 50, nobreak: true, readonly: true},
      { name: 'var_v_x', type: 'float', label: '&plusmn;', default: 1, min: 0.5, max: 2.0},
      { name: 'v_y', type: 'float', label: 'v<sub>y</sub>', default: 0, units: 'm/s' , readonly: true},
      { name: 'mean_w_x', type: 'float', label: 'w<sub>x</sub>', default: 0, nobreak: true, readonly: true},
      { name: 'var_w_x', type: 'float', label: '&plusmn;', default: 2, readonly: true  }
    ]
  }, {
    label: 'Mission Constraints', params: [
      { name: 'x_min', type: 'float', label: 'x<sub>min</sub>', default: -50, units: 'm', readonly: true },
      { name: 'x_max', type: 'float', label: 'x<sub>max</sub>', default: 50, units: 'm', readonly: true },
      { name: 'v_max', type: 'float', label: 'v<sub>max</sub>', default: 10, units: 'm/s', readonly: true },
      { name: 'F_max', type: 'float', label: 'F<sub>max</sub>', default: 300, units: 'N', readonly: true }
    ]
  }, {
    label: 'Monte Carlo', params: [
      { name: 'n_samp', type: 'int', label: 'Samples', default: 25000 }
    ]
  }
];

Simulation.prototype.getLabel = function (name) {
  this.metadata.forEach(function(group) {
    group.params.forEach(function(param, i) {
      if (name == param.name)
        return group.params[i].label;
    });
  });
};

Simulation.prototype.run = function (params) {
  var x   = params.x,    y      = params.y,
      v_x = params.v_x,  v_y    = params.v_y,
      r   = params.r,    r_p    = params.r_p,
      C_d = params.C_d,  C_dp   = params.C_dp,
      t_p = params.t_p,  t_open = params.t_open,
      w_x = params.w_x,
      m   = params.m;
  var dt = 0.02;
  var rho = 1.22, g = 9.8;
  var trajectory = [ ];
  var t = 0, p = 0, broken = 0; // parachute is closed at t=0
  while (y > 0) {
    if (params.save_trajectory === true)
      trajectory.push({t: t, x: x, y: y, v_x: v_x, v_y: v_y, p: p, b: broken});
    var v = Math.sqrt((v_x - w_x) * (v_x - w_x) + v_y * v_y);
    var F_d  = -0.5 * rho * C_d * (Math.PI * r * r) * v * v;
    var F_dp = -0.5 * rho * p * C_dp * (Math.PI * r_p * r_p) * v * v;
    if (Math.abs(F_dp) > params.F_max)
      broken = 1;
    var a_x = -Math.abs(F_d + F_dp) / m * (v_x - w_x) / v;
    var a_y = -g - Math.abs(F_d + F_dp) / m * v_y / v;
    var a = Math.sqrt(a_x * a_x + a_y * a_y);
    // dynamic step size
    dt = Math.min(0.5, 0.02 + 2 / (a + 1));
    if (Math.abs(t - t_p) < 1.1 * dt)
      dt = 0.04;
    // explicit euler time step
    v_x = v_x + a_x * dt;
    v_y = v_y + a_y * dt;
    x = x + v_x * dt;
    y = y + v_y * dt;
    t = t + dt;
    // update parachute open status
    if (t < t_p || broken) p = 0;
    else if (t > t_p + t_open) p = 1;
    else p = Math.min(1, (t - t_p) / t_open);
  }
  // lerp to y = 0
  x = x - y * v_x / v_y;
  trajectory.push({t: t, x: x, y: 0, v_x: v_x, v_y: v_y, p: p,  b: broken});
  return trajectory;
}

Simulation.prototype.getSamples = function(params) {
  var samples = {
    w_x: new Float64Array(params.n_samp),
    t_open: new Float64Array(params.n_samp),
    C_d: new Float64Array(params.n_samp),
    m: new Float64Array(params.n_samp),
    v_x: new Float64Array(params.n_samp)
  };
  // want log normal with mean and variance
  var mean = params.mean_t_open;
  var variance = params.var_t_open;
  // sample from Exp[N(mu, sigma2)]
  var lognormal_params = function(mean, variance) {
    var mu = Math.log(mean * mean / Math.sqrt(mean * mean + variance));
    var sigma2 = 2 * Math.log(Math.sqrt(mean * mean + variance) / mean);
    return {mu: mu, sigma2: sigma2};
  };
  var t_open_params = lognormal_params(params.mean_t_open, params.var_t_open);
  var m_params = lognormal_params(params.mean_m, params.var_m);
  for (var i = 0; i < params.n_samp; i++) {
    samples.w_x[i] = Simulation.getNormal() * Math.sqrt(params.var_w_x) + params.mean_w_x;
    samples.t_open[i] = Math.exp(Simulation.getNormal() * Math.sqrt(t_open_params.sigma2) + t_open_params.mu);
    samples.C_d[i] = params.C_d + (Math.random() * 2 - 1) * params.C_derr;
    samples.m[i] = Math.exp(Simulation.getNormal() * Math.sqrt(m_params.sigma2) + m_params.mu);
    samples.v_x[i] = Simulation.getNormal() * Math.sqrt(params.var_v_x) + params.mean_v_x;
  }
  return samples;
};

Simulation.getNormal = function() {
  var x, y, w;
  do {
    x = Math.random() * 2 - 1;
    y = Math.random() * 2 - 1;
    w = x * x + y * y;
  } while (w >= 1.0)
  return x * Math.sqrt(-2 * Math.log(w) / w);
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