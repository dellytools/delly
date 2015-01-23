// Copyright (c) 2014 Markus Hsi-Yang Fritz

var nBinsMax = 10000;
var nBinsMin = 2000;

var suave = function () {
  var my = {};

  var data = null;
  var outerWidth = 960;
  var margin = { top: 50, bottom: 20, left: 50, right: 50 };
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
                {n: nBinsMax},
                function (res) {
        my.data = res;
        $.getJSON('/calls/' + c, function (res) {
          my.data.calls = res;
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
      .tickFormat(d3.format('.0f'));

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
      .on("zoom", zoomed);

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
    var colors = {'INV': 'orange', 'DUP': 'DodgerBlue', 'DEL': '#333'};
    var arcWidthDefault = 2;
    var arcWidthBold = 4;
    var posFormat = d3.format(',d');

    drawArcs();

    function drawArcs() {
      console.log('arcs', xDepth.invert(0), xDepth.invert(my.arc.width));
      arcG.selectAll('path')
        .data(my.data.calls)
        .enter()
        .append('path')
        .filter(function (d) {
          return d['start'] >= xDepth.invert(0)
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
       .style('stroke', function (d) {return colors[d['type']];})
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
        return d['type']
          + ' ('
          + d['ct']
          + ') '
          + 'start: '
          + posFormat(d['start'])
          + ' end: '
          + posFormat(d['end']);
      });
    }
    function zoomed() {
      var sliceStart = Math.max(Math.floor(xDepth.invert(0)), 0);
      var sliceEnd = Math.min(Math.ceil(xDepth.invert(my.depth.width)),
                              my.data.chrom_len);

      // these are only valid/used if
      // sliceStart >= dataSeqStart && sliceEnd <= dataSeqEnd
      var binStart = Math.floor((sliceStart - dataSeqStart) / binSize);
      var binEnd = Math.floor((sliceEnd - dataSeqStart) / binSize);
      var nBins = binEnd - binStart + 1;

      console.log(sliceStart, sliceEnd, dataSeqStart, dataSeqEnd);
      console.log(binStart, binEnd, nBins);

      if (sliceStart >= dataSeqStart && sliceEnd <= dataSeqEnd
          && nBins >= nBinsMin && nBins >= nBinsMin) {
        depthG.select(".x.axis").call(xAxisDepth);
        redraw(depthCtx, binStart, binEnd);
      } else {
        console.log('GET new data');
        $.getJSON('/depth/' + my.sample1 + '/' + my.sample2 + '/' + my.chrom,
                  {start: sliceStart, end: sliceEnd, n: nBinsMax},
                  function (res) {
          my.data.ratios = res.ratios;
          binSize = Math.ceil((sliceEnd-sliceStart+1)  / my.data.ratios.length);
          dataSeqStart = sliceStart;
          dataSeqEnd = sliceEnd;
          depthG.select(".x.axis").call(xAxisDepth);
          redraw(depthCtx, 0, my.data.ratios.length-1);
        });
      }
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
      arcG.selectAll('path').remove();
      drawArcs();
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
  my.main = main;
  my.vis = vis;

  return my;
}();

$(suave.main.bind(null, '#vis'));
