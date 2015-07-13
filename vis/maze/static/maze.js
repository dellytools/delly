// 2014-2015 Markus Hsi-Yang Fritz

var maze = function () {
  var my = {};

  my.data = null;

  my.main = function (selector) {
    $('#footer-icon').click(function () {
      var revealing = parseInt($('footer').css('bottom')) < 0;
      $('footer').animate({
        bottom: revealing
                ? 0
                : -1 * parseInt($('footer').css('height'))
      }, 800, function () {
        $('#footer-icon').toggleClass('fa-chevron-circle-up');
        $('#footer-icon').toggleClass('fa-chevron-circle-down');
      });
    }); 

    // remove focus after being clicked
    $('.button-header').focus(function () {
      this.blur();
    }); 

    $('#match-info').click(function () {
      $('#matches-info').toggleClass('hide');
    });
    
    $('#checkbox-scale').click(function () {
      if ($('#checkbox-scale').prop('checked')) {
        $('#config-dim').prop('disabled', true);
      } else {
        $('#config-dim').prop('disabled', false);
      }
    });

    $('#drop-ref').on('dragenter', dropZoneEnter);
    $('#drop-ref').on('dragleave', dropZoneExit);
    $('#drop-ref').on('dragover', dropZoneOver);
    $('#drop-ref').on('drop', function (e) {
      $(this).removeClass('dropzone-hover');
      // TODO possible to removeClass "fa-*"?
      $('#ref-icon-chosen').removeClass('fa-times');
      $('#ref-icon-chosen').removeClass('fa-check');
      $('#ref-icon-chosen').addClass('fa-spinner').addClass('fa-pulse');
      getDroppedFasta(e, 'ref');
    });

    $('#drop-query').on('dragenter', dropZoneEnter);
    $('#drop-query').on('dragleave', dropZoneExit);
    $('#drop-query').on('dragover', dropZoneOver);
    $('#drop-query').on('drop', function (e) {
      $(this).removeClass('dropzone-hover');
      $('#query-icon-chosen').removeClass('fa-times');
      $('#query-icon-chosen').removeClass('fa-check');
      $('#query-icon-chosen').addClass('fa-spinner').addClass('fa-pulse')
;
      getDroppedFasta(e, 'query');
    });

    function dropZoneEnter(e) {
      e.stopPropagation();
      e.preventDefault();
      $(this).addClass('dropzone-hover');
    }

    function dropZoneExit(e) {
      e.stopPropagation();
      e.preventDefault();
      $(this).removeClass('dropzone-hover');
    }

    function dropZoneOver(e) {
      e.stopPropagation();
      e.preventDefault();
    } 

    function getDroppedFasta(e, type) {
      e.stopPropagation();
      e.preventDefault();

      var dt = e.dataTransfer || (e.originalEvent && e.originalEvent.dataTransfer);
      var files = e.target.files || (dt && dt.files);
      var f = files[0];

      readFasta(f, type, /\.gz$/.test(f.name));
    }

    function readFasta(f, type, isGzip) {
      var fReader = new FileReader();
      if (isGzip) {
        fReader.readAsArrayBuffer(f);
      } else {
        fReader.readAsText(f);
      }
      fReader.onload = function (e) {
        var fContent = e.target.result;
        if (isGzip) {
          fContent = pako.ungzip(fContent, {'to': 'string'});
        }

        var seqs = parseFastaString(fContent);

        if (seqs.length > 0) {
          my[type] = seqs;
          if (type === 'ref') {
            $('#ref-icon-chosen').removeClass('fa-spinner')
                                 .removeClass('fa-pulse')
                                 .addClass('fa-check');
            $('#ref-span-chosen').html(f.name);
          } else {
            $('#query-icon-chosen').removeClass('fa-spinner')
                                 .removeClass('fa-pulse')
                                 .addClass('fa-check');
            $('#query-span-chosen').html(f.name);
          }
        } else {
          if (type === 'ref') {
            $('#ref-icon-chosen').removeClass('fa-spinner')
                                 .removeClass('fa-pulse')
                                 .addClass('fa-times');
            $('#ref-span-chosen').empty();
          } else {
            $('#query-icon-chosen').removeClass('fa-spinner')
                                 .removeClass('fa-pulse')
                                 .addClass('fa-times');
            $('#query-span-chosen').empty();
          }
        }

        if (my.ref && my.query) {
          $('#visualize').prop('disabled', false);
        }
      }
    }

    $('#visualize').click(function () {
      $('#config-modal').modal('hide');
      $(selector).empty();
      $('#control-btn-left').addClass('hide');
      $('#control-btn-right').addClass('hide');
      $('.spinner').toggleClass('hide');

      var matches = $('#config-matches label.active').text().trim();
      var length = $('#config-length').val();

      $.post('/matches', {
          'matches': matches,
          'length': length,
          'ref': JSON.stringify(my.ref[0]),
          'query': JSON.stringify(my.query)
        }, function (res) {
          $('.spinner').toggleClass('hide');
          my.data = res;
          my.vis(selector, 0);
        }, 'json'
      );
    });
  };

  my.vis = function (selector, dataIdx) {
    var data = my.data[dataIdx]
    var l1 = my.ref[0].seq.length;
    var l2 = my.query[dataIdx].seq.length;

    $(selector).empty();

    $('#control-btn-left').removeClass('hide');
    $('#control-btn-right').removeClass('hide');
    $('#control-btn-breakpoints').off('click');

    $('#control-btn-breakpoints:not(disabled)').click(function () { // currently never disabled
        console.log('open new window for ' + dataIdx)
        var wnd = window.open("breakpoints");
        wnd.transferData = { data: data, query: my.query[dataIdx], ref: my.ref[0]}; // Todo(meiers): Change ref[0] in the future
      });

    if (dataIdx > 0) {
      $('#control-btn-left').removeClass('disabled');
      $('#control-btn-left:not(disabled)').one('click', function () {
        my.vis(selector, dataIdx-1);
      });
    } else {
      $('#control-btn-left').addClass('disabled');
    }

    if (dataIdx < my.data.length - 1) {
      $('#control-btn-right').removeClass('disabled');
      $('#control-btn-right:not(disabled)').one('click', function () {
        my.vis(selector, dataIdx+1);
      });
    } else {
      $('#control-btn-right').addClass('disabled');
    }

    if ($('#checkbox-scale').prop('checked')) {
      my.outerWidth = Math.min($(window).width(),
                               $(window).height()) * 0.8;
    } else {
      my.outerWidth = +$('#config-dim').val();
    }

    my.outerHeight = my.outerWidth;
    my.margin = { top: 25, bottom: 15, left: 50, right: 50 };
    var innerWidth = my.outerWidth - my.margin.left - my.margin.right;
    my.innerWidth = innerWidth;
    var innerHeight = my.outerHeight - my.margin.top - my.margin.bottom;
    my.innerHeight = innerHeight;

    var x = d3.scale.linear()
      .domain([1, l1])
      .range([0, my.innerWidth]);

    var y = d3.scale.linear()
      .domain([1, l2])
      .range([0, my.innerHeight]);

    var xAxisB = d3.svg.axis()
      .scale(x)
      .orient('bottom')
      .tickSize(-my.innerHeight);

    var xAxisT = d3.svg.axis()
      .scale(x)
      .orient('top')
      .tickSize(0);

    var yAxisR = d3.svg.axis()
      .scale(y)
      .orient('right')
      .tickSize(-my.innerWidth);

    var yAxisL = d3.svg.axis()
      .scale(y)
      .orient('left')
      .tickSize(-my.innerWidth);
  
    var zoom = d3.behavior.zoom()
      .x(x)
      .y(y)
      .scaleExtent([1, 100])
      .on('zoom', zoomed);

    var svg = d3.select(selector).append('svg')
      .attr('width', my.outerWidth)
      .attr('height', my.outerHeight);

    var g = svg.append('g')
      .attr('transform', 'translate(' + my.margin.left + ', ' + my.margin.top + ')');

    g.call(zoom);

    g.append('rect')
      .attr('class', 'overlay')
      .attr('width', my.innerWidth)
      .attr('height', my.innerHeight);

    g.append('g')
      .attr('class', 'x axis b')
      .attr('transform', 'translate(0, ' + my.innerHeight + ')')
      .call(xAxisB);

    g.append('g')
      .attr('class', 'x axis t')
      .call(xAxisT);

    g.append('g')
      .attr('class', 'y axis r')
      .attr('transform', 'translate(' + my.innerWidth + ', 0)')
      .call(yAxisR);

    g.append('g')
      .attr('class', 'y axis l')
      .call(yAxisL);

    g.append('defs')
      .append('svg:clipPath')
      .attr("id", "clipChartArea")
      .append('rect')
      .attr('x', 0)
      .attr('y', 0)
      .attr('width', my.innerWidth)
      .attr('height', my.innerHeight);

    var chartArea = g.append('g')
      .attr('clip-path', 'url(#clipChartArea)');

    chartArea.selectAll('line.matches')
      .data(data.fwd.concat(data.rev))
      .enter()
      .append('line')
      .attr('class', 'geom matches')
      .attr('x1', function (d) { return x(d[0]); })
      .attr('x2', function (d) { return x(d[1]); })
      .attr('y1', function (d) { return y(d[2]); })
      .attr('y2', function (d) { return y(d[3]); })
      // TODO should do this via CSS classes:
      .style('stroke', function (d) {
        return d[2] < d[3] ? 'black' : 'red'
      });

      function zoomed() {
        chartArea.selectAll('.geom')
          .data(data.fwd.concat(data.rev))
          .attr('x1', function (d) { return x(d[0]); })
          .attr('x2', function (d) { return x(d[1]); })
          .attr('y1', function (d) { return y(d[2]); })
          .attr('y2', function (d) { return y(d[3]); });
        g.select(".x.axis.b").call(xAxisB);
        g.select(".x.axis.t").call(xAxisT);
        g.select(".y.axis.r").call(yAxisR);
        g.select(".y.axis.l").call(yAxisL);
      }
  };

  return my;
}();

$(maze.main.bind(null, '#vis'));

function LineReader(str) {
  this.str = str;
  this.i = 0;
  this.next = nextLine;
}

function nextLine() {
  var l = '';
  for (; this.i < this.str.length; this.i += 1) {
    if (this.str[this.i] === '\n') {
      this.i += 1;
      break;
    }
    l += this.str[this.i];
  }
  return l;
}

function parseFastaString(s) {
  var seqs = [];
  var name = null;
  var seq = null;
  var lr;
  var l;

  lr = new LineReader(s);
  while (l = lr.next()) {
    if (l[0] === '>') {
      if (name) {
        seqs.push({'name': name, 'seq': seq});
      }
      seq = '';
      name = />(\w+)/.exec(l)[1];
    } else {
      seq += l;
    }
  }

  if (name) {
    seqs.push({'name': name, 'seq': seq});
  }

  return seqs;
}
