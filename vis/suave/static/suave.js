// Copyright (c) 2014 Markus Hsi-Yang Fritz

var nBinsMax = 20000;
var nBinsMin = Math.ceil(nBinsMax / 10);

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
      $.getJSON('/data/' + s1 + '/' + s2 + '/' + c,
                {n: nBinsMax},
                function (res) {
        my.data = res;
        my.vis(selector);
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
    var viewStart = 1;
    var viewEnd = my.data.chrom_len;

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

    function zoomed() {
      var sliceStart = Math.max(Math.floor(xDepth.invert(0)), 0);
      var binStart = Math.floor(sliceStart / binSize);
      var sliceEnd = Math.min(Math.ceil(xDepth.invert(my.depth.width)),
                              my.data.chrom_len);
      var binEnd = Math.floor(sliceEnd / binSize);
      var i;

      console.log(sliceStart, sliceEnd);
      console.log(binStart, binEnd, binEnd - binStart + 1);

      if (binEnd - binStart + 1 > nBinsMin) {
        depthG.select(".x.axis").call(xAxisDepth);
        redrawCanvas(depthCtx, binStart, binEnd);
      } else {
        $.getJSON('/data/' + my.sample1 + '/' + my.sample2 + '/' + my.chrom,
                  {start: sliceStart, end: sliceEnd, n: nBinsMax},
                  function (res) {
          my.data = res;
          binSize = Math.ceil((sliceEnd-sliceStart+1)  / my.data.ratios.length);
          depthG.select(".x.axis").call(xAxisDepth);
          redrawCanvas(depthCtx, 0, my.data.ratios.length-1);
        });
      }
    }

    function redrawCanvas(ctx, start, end) {
      var i;
      ctx.clearRect(0, 0, my.depth.width, my.depth.height);
      for (i = start; i <= end; i += 1) {
        canvasPoint(ctx, xDepth(i*binSize), yDepth(my.data.ratios[i]));
      }
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
