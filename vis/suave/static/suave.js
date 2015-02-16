// Copyright (c) 2014 Markus Hsi-Yang Fritz

// TODO: we need to benchmark this
// also: at some point we could dynamically set this based on 
//       screen properties
var N_BINS_DESIRED = 10000;

var suave = function () {
  var my = {};

  var data = null;
  var outerWidth = 960;
  var margin = { top: 10, bottom: 20, left: 50, right: 50 };
  var innerWidth = outerWidth - margin.left - margin.right;
  var arc = {
    width: innerWidth,
    height: 150,
    offX: margin.left,
    offY: margin.top
  };
  var depth = {
    width: innerWidth,
    height: 250,
    offX: margin.left,
    offY: arc.offY + arc.height
  };
  var brush = {
     width: innerWidth,
     height: 75,
     offX: margin.left,
     offY: depth.offY + depth.height
  };
  var innerHeight = arc.height + depth.height + brush.height;
  var outerHeight = innerHeight + margin.top + margin.bottom;

  // TODO this needs to be dynamic based on the actual VCF
  // also: use colorbrewer palette
  var svColors = {'INV': 'orange', 'DUP': '#666', 'DEL': 'DodgerBlue'};

  var main = function (selector) {
    $('.sample_select').click(function () {
      var s = {
        'sample1': $('#sample1_selected'),
        'sample2': $('#sample2_selected')
      };
      var selected = $(this).data("sample");
      var other = selected === 'sample1' ? 'sample2' : 'sample1';
      var s1, s2;

      s[selected].html($(this).html());

      // FIXME
      if (s[other].html() !== '-select-') {
        s1 = $('#sample1_selected').html();
        s2 = $('#sample2_selected').html();
        $.getJSON('/chroms/' + s1 + '/' + s2, function (res) {
          var select_html = '';
          $.each(res, function (i, v) {
            select_html += '<li><a href="#" class="chromosome_select">' +
                             v + '</a>' + '</li>';
          });
          $('#chromosome_dropdown>ul').html(select_html);
          $('.chromosome_select').click(function () {
            $('#chromosome_selected').html($(this).html());
          });
        });
      }
    });

    $('#runButton').click(function () {
      var s1, s2, c;
      s1 = $('#sample1_selected').html();
      s2 = $('#sample2_selected').html();
      c = $('#chromosome_selected').html();
      // FIXME
      if (s1 === '-select-' || s2 === '-select-' || c === '-select-') {
        if ($('.alert-danger').length === 0) {
          $('#selectedArea').after('<div class="alert alert-danger">' +
                                     'Please specify samples and chromosome!' +
                                   '</div>');
        }
        return;
      }
      $('.alert-danger').remove();
      my.sample1 = s1;
      my.sample2 = s2;
      my.chrom = c;
      $.getJSON('/depth/' + s1 + '/' + s2 + '/' + c,
                {n: N_BINS_DESIRED},
                function (res) {
        my.data = res;
        $.getJSON('/calls/' + c, function (res) {
          my.data.calls = res;
          var svTypes = d3.set(_.pluck(suave.data.calls, 'type'))
                          .values()
                          .sort();
          $('#svTypes').remove();
          $('#vis-wrap').prepend('<div id="svTypes"></div>');
          $.each(svTypes, function (idx, val) {
            var html = '<input type="checkbox" checked value="'
                       + val + '"> <span style="color:' 
                       + svColors[val] + '">' + val + ' </span>';
            $('#svTypes').append(html);
          });

          $('#jumpToSlice').removeClass('hide');
          $('#jumpToPad').val('0');

          $('#jumpToFeature').removeClass('hide');
          $('#jumpToPadFeature').val('100');

          var bloodhound = new Bloodhound({
            datumTokenizer: Bloodhound.tokenizers.obj.whitespace('value'),
            queryTokenizer: Bloodhound.tokenizers.whitespace,
            limit: 100,
            local: $.map(my.data.calls, function(call) { return { value: call.id }; })
          });

          bloodhound.initialize();

          $('#bloodhound .typeahead').typeahead(null, {
            name: 'bloodhound',
            displayKey: 'value',
            source: bloodhound.ttAdapter()
          });

          my.vis(selector);
        });
      });
    });
  };

  var vis = function (selector) {
    var i;
    var circleRadians = 2 * Math.PI;

    $(selector).empty();

    var frameSvg = d3.select(selector).append('svg')
      .attr('width', my.outerWidth)
      .attr('height', my.outerHeight);

    var arcG = frameSvg.append('g')
      .attr('transform', 'translate(' + my.arc.offX
            + ', ' + my.arc.offY + ')');
    var depthG = frameSvg.append('g')
      .attr('transform', 'translate(' + my.depth.offX
            + ', ' + my.depth.offY + ')');
    var brushG = frameSvg.append('g')
      .attr('transform', 'translate(' + my.brush.offX
            + ', ' + my.brush.offY + ')');

    // ------------ Read Depth ---------------------

    var depthCanvas = d3.select(selector).append('canvas')
      .attr('width', my.depth.width)
      .attr('height', my.depth.height)
      .style('left', String(my.depth.offX) + 'px')
      .style('top', String(my.depth.offY) + 'px');
    var depthCtx = depthCanvas.node().getContext('2d');

    var xDepth = d3.scale.linear()
      .domain([1, my.data.chrom_len])
      .range([0, my.depth.width]);

    var yDepth = d3.scale.linear()
      .domain(d3.extent(my.data.ratios))
      .range([my.depth.height, 0]);

    var xAxisDepth = d3.svg.axis()
      .scale(xDepth)
      .orient('bottom')
      .tickSize(-my.depth.height);

    depthG.append('g')
      .attr('class', 'x axis')
      .attr('transform', 'translate(0, ' + my.depth.height + ')')
      .call(xAxisDepth);

    var yAxisDepth = d3.svg.axis()
      .scale(yDepth)
      .orient('left')
      .tickSize(-my.depth.width)
      .tickFormat(d3.format('.1f'));

    depthG.append('g')
      .attr('class', 'y axis')
      .call(yAxisDepth);

    var binSize = Math.ceil(my.data.chrom_len / my.data.ratios.length);
    var dataSeqStart = 1;
    var dataSeqEnd = my.data.chrom_len;

    $.each(my.data.ratios, function (idx, val) {
      canvasPoint(depthCtx, xDepth(idx*binSize), yDepth(val));
    });

    var zoom = d3.behavior.zoom()
      .x(xDepth)
      .scaleExtent([1, Infinity])
      .on("zoom", rescale.bind(null, 'zoom'));

    var zoomLocked = false;

    depthCanvas.call(zoom);

    function canvasPoint(c, x, y) {
      c.beginPath();
      c.arc(x, y, 1, 0, circleRadians, true);
      c.fill();
    };

    // ------------ Arcs ---------------------
  
    var arc = d3.svg.line()
      .x(function(d) { return xDepth(d.x); })
      .y(function(d) { return d.y; })
      .interpolate('cardinal');
    var arcWidthDefault = 2;
    var arcWidthBold = 4;
    var posFormat = d3.format(',d');

    $('input').on('click', drawArcs);
    drawArcs();

    function drawArcs() {
      arcG.selectAll('path').remove();

      var checked = d3.map();
      $.each($('input'), function (idx, val) {
        checked.set(val.value, val.checked);
      });

      arcG.selectAll('path')
        .data(my.data.calls)
        .enter()
        .append('path')
        .filter(function (d) {
          return checked.get(d['type'])
              && d['start'] >= xDepth.invert(0)
              && d['end'] <= xDepth.invert(my.arc.width);
        })
        .attr('d', function (d) {
          return arc([
            {x: d['start'], y: my.arc.height},
            {x: (d['start'] + d['end']) / 2, y: 2},
            {x: d['end'], y: my.arc.height}
         ]);
       })
       .style('fill', 'none')
       .style('stroke', function (d) {return svColors[d['type']];})
       .style('stroke-width', arcWidthDefault)
       .on('mouseover', function () {
        arcG.selectAll('path')
          .style('opacity', '0.1');
        d3.select(this)
         .style('opacity', '1')
         .style('stroke-width', arcWidthBold);
      })
      .on('mouseout', function () {
        arcG.selectAll('path')
          .style('opacity', '1');
        d3.select(this)
          .style('stroke-width', arcWidthDefault);
      })
      .append('title').text(function (d) {
        return d['id']
          + '\x0A'
          + d['type']
          + ' ('
          + d['ct']
          + ')\x0A'
          + 'start: '
          + posFormat(d['start'])
          + ' end: '
          + posFormat(d['end']);
      });
    }

    // ------------ Brush ---------------------

    var brushCanvas = d3.select(selector).append('canvas')
      .attr('width', my.brush.width)
      .attr('height', my.brush.height)
      .style('left', String(my.brush.offX) + 'px')
      .style('top', String(my.brush.offY) + 'px')
      .style('z-index', '-10');
    var brushCtx = brushCanvas.node().getContext('2d');

    var xBrush = d3.scale.linear()
      .domain([1, my.data.chrom_len])
      .range([0, my.brush.width]);

    var xAxisBrush = d3.svg.axis()
      .scale(xBrush)
      .orient('bottom');

    // FIXME: hack to avoid overflowing depth <g>
    var yBrush = d3.scale.linear()
      .domain(d3.extent(my.data.ratios))
      .range([my.brush.height, 20]);
      //.range([my.brush.height, 0]);

    brushG.append('g')
      .attr('class', 'x axis')
      .attr('transform', 'translate(0, ' + my.brush.height + ')')
      .call(xAxisBrush);

    var brush = d3.svg.brush()
      .x(xBrush)
      .on('brush', rescale.bind(null, 'brush'));

    // FIXME: hack to avoid overflowing depth <g>
    var brusher = brushG.append('g')
      .attr('class', 'x brush')
      .call(brush);

    brusher.selectAll('rect')
      //.attr('y', -6)
      .attr('y', -6 + 20)
      //.attr('height', my.brush.height + 7);
      .attr('height', my.brush.height + 7 - 20);

    $.each(my.data.ratios, function (idx, val) {
      canvasPoint(brushCtx, xBrush(idx*binSize), yBrush(val));
    });

    $('#jumpToSliceSubmit').click(function () {
      var s = parseInt($('#jumpToStart').val());
      var e = parseInt($('#jumpToEnd').val());
      var p = parseInt($('#jumpToPad').val());

      if (isNaN(s) || isNaN(e) || isNaN(p)) {
        console.log('error: jumpToSlice NaN', s, e, p);
        return;
      }

      s -= p;
      e += p;

      s = s < 1 ? 1 : s;
      e = e > my.data.chrom_len ? my.data.chrom_len : e;

      if (s > my.data.chrom_len || e < 1 || s > e) {
        console.log('error: jumpToSlice coords', s, e, p);
        return;
      }

      rescale('jump', s, e);
    });

    $('#jumpToFeatureSubmit').click(function () {
      var featID = $('#jumpToFeatureInput').val();
      var p = parseInt($('#jumpToPadFeature').val());

      if (isNaN(p)) {
        console.log('error: jumpToFeature NaN', p);
        return;
      }

      // FIXME should index this...
      $.each(my.data.calls, function (idx, val) {
        if (val.id === featID) {
          console.log(val);
          var s = val.start - p;
          var e = val.end + p;
          s = s < 1 ? 1 : s;
          e = e > my.data.chrom_len ? my.data.chrom_len : e;
          rescale('jump', s, e);
          return false;
        }
      });
    });

    function rescale(control, start, end) {
      if (zoomLocked) {
        return;
      }

      zoomLocked = true;

      if (control === 'zoom') {
        brusher.call(brush.extent(xDepth.domain()))
          .call(brush.event);
      } else if (control === 'brush') {
        xDepth.domain(brush.empty() ? xBrush.domain() : brush.extent()); 
        zoom.x(xDepth);
      } else if (control === 'jump') {
        xDepth.domain([start, end]);
        zoom.x(xDepth);
        brusher.call(brush.extent(xDepth.domain()))
          .call(brush.event);
      }

      var sliceStart = Math.max(Math.floor(xDepth.invert(0)), 1);
      var sliceEnd = Math.min(Math.ceil(xDepth.invert(my.depth.width)),
                              my.data.chrom_len);

      console.log(sliceStart, sliceEnd, dataSeqStart, dataSeqEnd);

      console.log('GET new data');
      $.getJSON('/depth/' + my.sample1 + '/' + my.sample2 + '/' + my.chrom,
                {start: sliceStart, end: sliceEnd, n: N_BINS_DESIRED},
                function (res) {
        my.data.ratios = res.ratios;
        binSize = Math.ceil((sliceEnd-sliceStart+1)  / my.data.ratios.length);
        dataSeqStart = sliceStart;
        dataSeqEnd = sliceEnd;
        depthG.select(".x.axis").call(xAxisDepth);
        redraw(depthCtx, 0, my.data.ratios.length-1);
      });
    }

    function redrawCanvas(ctx, start, end) {
      var i;
      ctx.clearRect(0, 0, my.depth.width, my.depth.height);
      for (i = start; i <= end; i += 1) {
        canvasPoint(ctx, xDepth(i*binSize+dataSeqStart), yDepth(my.data.ratios[i]));
      }
    }

    function redraw(ctx, start, end) {
      redrawCanvas(ctx, start, end);
      drawArcs();
      zoomLocked = false;
    }
  };

  my.data = data;
  my.outerWidth = outerWidth;
  my.margin = margin;
  my.innerWidth = innerWidth;
  my.arc = arc;
  my.depth = depth;
  my.brush = brush;
  my.innerHeight = innerHeight;
  my.outerHeight = outerHeight;
  my.svColors = svColors;
  my.main = main;
  my.vis = vis;

  return my;
}();

$(suave.main.bind(null, '#vis'));
