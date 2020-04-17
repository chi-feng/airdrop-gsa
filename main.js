'use strict';

var sim_worker = new Worker('simulation-worker.js');

var plots = { };

var hist = function(data, min, max, n_bins, value, filter) {
  var width = (max - min) / (n_bins - 1);
  var counts = [];
  var bins = [];
  for (var i = 0; i < n_bins; i++) {
    bins.push(i * width + min);
    counts.push(0);
  }
  for (var i = 0; i < data.params.n_samp; i++) {
    var bindex = Math.floor((value(data, i) - min) / width);
    if (bindex >= 0 && bindex < n_bins && filter(data, i))
      counts[bindex]++;
  }
  return {bins: bins, counts: counts};
};

var run_sa_handler = function(event) {
  var keys = event.data.keys;
  var S = event.data.main_effect_sensitivites;
  var ST = event.data.total_effect_sensitivites;

  function pie_data(keys, sens) {
    var values = [], labels=[];
    keys.forEach(function(key, j) {
      values.push(sens[j]);
      labels.push(key);
    });
    var sum = 0;
    for (var j = 0; j < keys.length; j++)
      sum += sens[j];
    values.push(1 - sum);
    labels.push('H/O terms');
    return {values: values, labels: labels, type: 'pie', sort: false,
      marker: {colors:[
        'rgb(102, 194, 165)',
        'rgb(252, 141, 98)',
        'rgb(141, 160, 203)',
        'rgb(231, 138, 195)',
        'rgb(166, 216, 84)',
        'rgb(240, 240, 240)',
        'rgb(229, 196, 148)',
        'rgb(255, 217, 47)',
      ]}
    };
  };

  function bar_data(keys, sens) {
    var x = [], y = [];
    keys.forEach(function(key, j) {
      y.push(sens[j]);
      x.push(key);
    });
    return {x: x, y: y, type: 'bar'};
  };

  var x_main = pie_data(keys, S.x);
  var x_total = bar_data(keys, ST.x);
  var v_main = pie_data(keys, S.v);
  var v_total = bar_data(keys, ST.v);

  x_main.domain = { x: [0, .48] };
  x_main.hole = 0.4;
  v_main.domain = { x: [.52, 1] };
  v_main.hole = 0.4;

  x_total.name = 'x';
  v_total.name = 'v';

  var data = [x_main, v_main];

  Plotly.newPlot(plots.xSensPlot, data, {
    margin: {l: 20, r: 20, b: 20, t: 50 }, title: 'Main effect sensitivities',
    titlefont: { size:16, family:'Helvetica, Arial' },
    annotations: [{
      font: { size: 20 },
      showarrow: false,
      text: 'x(t<sub>f</sub>)',
      x: 0.17, y: 0.5
    }, {
      font: { size: 20 },
      showarrow: false,
      text: 'v(t<sub>f</sub>)',
      x: 0.82, y: 0.5
    }]
  });

  Plotly.newPlot(plots.xTotalPlot, [x_total, v_total], {
    margin: {l: 40, r: 20, b: 30, t: 50 }, title: 'Total effect sensitivities',
    titlefont: { size:16, family:'Helvetica, Arial' }
  });

  reset_button();
};

var reset_button = function() {
  var btn = document.getElementById('btn_run');
  btn.disabled = false;
  btn.innerHTML = 'Run Monte Carlo';
}

var run_mc_handler = function(event) {
  var data = event.data;
  var params = data.params;
  // run SA
  var btn = document.getElementById('btn_run');
  btn.innerHTML = 'Running Sensitivity Analysis';
  params.action = 'run_sa';
  sim_worker.postMessage(params);
  var html = '';
  var success = 0;
  var outside = 0;
  var torn = 0;
  var destroyed = 0;
  var unopened = 0;
  for (var i = 0; i < params.n_samp; i++) {
    if (data.v[i] > params.v_max)
      destroyed++;
    if (data.x[i] >= params.x_min && data.x[i] <= params.x_max) {
      if (data.v[i] <= params.v_max)
        success++;
    } else {
      if (data.v[i] <= params.v_max)
        outside++;
    }
    if (data.p[i] < 1)
      unopened++;
    if (data.b[i] == 1)
      torn++;
  }
  html += '<strong>Successful landing:</strong> ' + parseFloat(100 * success / params.n_samp).toFixed(1) + '%<br />';
  html += '<strong>Destroyed on impact:</strong> ' + parseFloat(100 * destroyed / params.n_samp).toFixed(1) + '%<br />';
  html += '<strong>Landed outside:</strong> ' + parseFloat(100 * outside / params.n_samp).toFixed(1) + '%<br />';
  // html += '<strong>Parachute torn off:</strong> ' + parseFloat(100 * torn / params.n_samp).toFixed(1) + '%<br />';
  // html += '<strong>Parachute not fully open:</strong> ' + parseFloat(100 * unopened / params.n_samp).toFixed(1) + '%';

  document.getElementById('summary').innerHTML = html;

  var layout = {
    yaxis: {
        showgrid: true,
        zeroline: true,
        dtick: 5,
        gridcolor: 'rgba(0,0,0,0.1)',
        gridwidth: 1,
        zerolinecolor: 'rgba(0,0,0,0.25)',
        zerolinewidth: 1
    },
    margin: {l: 20, r: 20, b: 20, t: 20 },
    showlegend: false
  };

  var traces = [];
  data.trajectories.forEach(function(traj, i) {
    var before    = { x: [ ], y: [ ], mode: 'lines', type: 'scatter', line: { width: 1, color: 'rgba(240,120,0,0.5)' } };
    var after    = { x: [ ], y: [ ], mode: 'lines', type: 'scatter', line: { width: 1, color: 'rgba(30,100,200,0.5)' } };
    var open     = { x: [ ], y: [ ], mode: 'markers', type: 'scatter', marker: { symbol:'cross', color: 'rgba(30,100,200,1)', size: 4 } };
    var full     = { x: [ ], y: [ ], mode: 'markers', type: 'scatter', marker: { symbol:'circle', color: 'rgba(30,100,200,1)', size: 3 } };
    var torn     = { x: [ ], y: [ ], mode: 'markers', type: 'scatter', marker: { symbol:'x', color: 'rgba(200,40,30,0.9)', size: 6 } };
    var success  = { x: [ ], y: [ ], mode: 'markers', type: 'scatter', marker: { symbol:'circle', color: 'rgba(30,200,40,0.5)', size: 6 } };
    var fail     = { x: [ ], y: [ ], mode: 'markers', type: 'scatter', marker: { symbol:'x', color: 'rgba(200,40,30,0.9)', size: 6 } };
    before.x = traj.before.x;
    before.y = traj.before.y;
    after.x = traj.after.x;
    after.y = traj.after.y;
    traces.push(before);
    traces.push(after);
    if (traj.success) {
      success.x.push(traj.land.x);
      success.y.push(traj.land.y);
      traces.push(success);
    } else {
      fail.x.push(traj.land.x);
      fail.y.push(traj.land.y);
      traces.push(fail);
    }
    if (traj.torn.y > -1) {
      torn.x.push(traj.torn.x);
      torn.y.push(traj.torn.y);
      traces.push(torn);
    }
    if (traj.full.y > -1) {
      full.x.push(traj.full.x);
      full.y.push(traj.full.y);
      traces.push(full);
    }
    if (traj.open.y > -1) {
      open.x.push(traj.open.x);
      open.y.push(traj.open.y);
      traces.push(open);
    }
  });

  var hist_min = params.x_min - (params.x_max - params.x_min) / 2;
  var hist_max = params.x_max + (params.x_max - params.x_min) / 2;
  var n_bins = 101;

  Plotly.newPlot(plots.trajPlot, traces, {

    // xaxis: {range: [hist_min, hist_max]},
    // yaxis: {range: [-10, hist_max - hist_min]},
    margin: {l: 30, r: 20, b: 20, t: 50 },
    showlegend: false,
    title: 'Simulated trajectories (showing first 25)',
    titlefont: { size:16, family:'Helvetica, Arial' }
  });


  var success_hist = hist(data, hist_min, hist_max, n_bins,
    function(data, i) {
      return data.x[i];
    },
    function(data, i) {
      return data.x[i] > params.x_min && data.x[i] < params.x_max && data.v[i] <= params.v_max;
    }
  );

  var destroy_hist = hist(data, hist_min, hist_max, n_bins,
    function(data, i) {
      return data.x[i];
    },
    function(data, i) {
      return data.v[i] > params.v_max;
    }
  );

  var outside_hist = hist(data, hist_min, hist_max, n_bins,
    function(data, i) {
      return data.x[i];
    },
    function(data, i) {
      return (data.x[i] < params.x_min || data.x[i] > params.x_max) && data.v[i] <= params.v_max;
    }
  );

  var success_trace = {
    x: success_hist.bins,
    y: success_hist.counts,
    name: 'Success',
    fill: 'tozeroy',
    fillcolor: 'rgba(30,90,190,0.5)',
    line: {shape: 'hv'},
    type: 'scatter',
    mode: 'none'
  };

  for (var i = 0; i < destroy_hist.counts.length; i++) {
    // destroy_hist.counts[i] += (success_hist.counts[i] + outside_hist.counts[i]);
  }

  var destroy_trace = {
    x: destroy_hist.bins,
    y: destroy_hist.counts,
    name: 'Destroyed',
    fill: 'tonexty',
    fillcolor: 'rgba(250,100,0,1.0)',
    line: {shape: 'hv'},
    type: 'scatter',
    mode: 'none'
  };

  var outside_trace = {
    x: outside_hist.bins,
    y: outside_hist.counts,
    name: 'Outside',
    fill: 'tozeroy',
    fillcolor: 'rgba(160,160,160,0.5)',
    line: {shape: 'hv'},
    type: 'scatter',
    mode: 'none'
  };

  Plotly.newPlot(plots.xPlot, [destroy_trace, success_trace, outside_trace, ], {
    margin: {l: 40, r: 20, b: 20, t: 50 },
    title: 'Simulated outcomes',
    titlefont: { size:16, family:'Helvetica, Arial' }
  });

};

sim_worker.onmessage = function(event) {
  if (event.data.action == 'run_mc') {
    run_mc_handler(event);
  }
  if (event.data.action == 'run_sa') {
    run_sa_handler(event);
  }
  if (event.data.action == 'error') {
    var errors = event.data.errors;
    window.alert(errors.join('\n'));
    reset_button();
  }
};

function getParams(id, metadata) {
  var params = { };
  metadata.forEach(function(group) {
    group.params.forEach(function(p) {
      var input = document.querySelector('input[name="' + p.name + '"]');
      params[p.name] = input.value;
    });
  });
  return params;
}

function getFluidPlot(id) {
  var d3 = Plotly.d3;
  var gd3 = d3.select('#' + id)
    .style({
      width: '100%', 'margin-left': '0%',
      height: '100%', 'margin-top': '0%',
    });
  return gd3.node();
}

function createForm(id, metadata) {
  var html = '';
  metadata.forEach(function(group) {
    html += '<fieldset><legend>' + group.label + '</legend>';
    group.params.forEach(function(p, i) {
      if (i > 0 && !group.params[i-1].nobreak)
        html += '<div class="param">';
      html += '<label>' + p.label + '</label>';
      var readonly = p.readonly ? 'readonly class="readonly"' : '';
      html += '<input type="text" name="' + p.name + '" value="' + p.default + '" ' + readonly + ' />';
      if (p.units != undefined)
        html += '<span class="units">' + p.units + '</span>';
      if (!p.nobreak)
        html += '</div>';
    });
    html += '</fieldset>';
  });
  html += '<fieldset><button id="btn_run">Run Monte Carlo</button></fieldset>';
  document.getElementById(id).innerHTML = html;
  var btn = document.getElementById('btn_run');
  btn.addEventListener('click', function(event) {
    event.preventDefault();
    btn.disabled = true;
    btn.innerHTML = 'Running Monte Carlo';
    var params = getParams(id, metadata);
    params.action = 'run_mc';
    sim_worker.postMessage(params);
  });
}

function resize(fluid) {

  var width = window.innerWidth;
  var height = window.innerHeight;

  var formWidth = document.getElementById('parameters').offsetWidth;
  var colWidth = Math.floor((width - formWidth) * 0.5) - 5;
  document.getElementById('col-traj').style.width  = colWidth + 'px';
  document.getElementById('col-traj').style.height = colWidth + 'px';
  document.getElementById('col-hist').style.width  = colWidth + 'px';
  document.getElementById('col-hist').style.height = colWidth + 'px';

  document.getElementById('plt-traj-container').style.width  = colWidth + 'px';
  document.getElementById('plt-traj-container').style.height = (colWidth * 0.8) + 'px';
  document.getElementById('plt-x-container').style.width  = colWidth + 'px';
  document.getElementById('plt-x-container').style.height = (colWidth * 0.65) + 'px';

  document.getElementById('plt-xsense-container').style.width  = colWidth + 'px';
  document.getElementById('plt-xsense-container').style.height = colWidth / 2 + 'px';
  document.getElementById('plt-xtotal-container').style.width  = colWidth + 'px';
  document.getElementById('plt-xtotal-container').style.height = colWidth / 2 + 'px';
  if (fluid) {
    Plotly.Plots.resize(plots.trajPlot);
    Plotly.Plots.resize(plots.xPlot);
    Plotly.Plots.resize(plots.xSensPlot);
    Plotly.Plots.resize(plots.xTotalPlot);
  }

}

function init() {
  var sim = new Simulation();
  createForm('parameters', sim.metadata);
  plots.trajPlot = getFluidPlot('plt-traj');
  plots.xPlot = getFluidPlot('plt-x');
  plots.xSensPlot = getFluidPlot('plt-xsense');
  plots.xTotalPlot = getFluidPlot('plt-xtotal');
  resize(false);
  window.onresize = function() { resize(true); };
}