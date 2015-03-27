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
        if (revealing) {
          $('#footer-icon').removeClass('fa-chevron-circle-up');
          $('#footer-icon').addClass('fa-chevron-circle-down');
        } else {
          $('#footer-icon').removeClass('fa-chevron-circle-down');
          $('#footer-icon').addClass('fa-chevron-circle-up');
        }
      });
    }); 

    // remove focus after being clicked
    $('.button-header').focus(function () {
      this.blur();
    }); 

    $('#configModal .modal-body .fa').click(function () {
      $('#matches-info').toggleClass('hide');
    });
    
    $('#visualize').click(function () {
      $('#configModal').modal('hide');
      $(selector).empty();
      $('#control-btn-left').addClass('hide');
      $('#control-btn-right').addClass('hide');
      $('.spinner').toggleClass('hide');

      var matches =$('#config-matches label.active').text().trim();
      var length = $('#config-length').val();

      $.getJSON('/matches',
        {'matches': matches, 'length': length},
        function (res) {
          $('.spinner').toggleClass('hide');
          my.data = res;
          my.vis(selector, 0);
        }
      );
    });
  };

  my.vis = function (selector, dataIdx) {
    $(selector).empty();

    var data = my.data[dataIdx]
    var l1 = data.rseq.length;
    var l2 = data.qseq.length;

    $('#control-btn-left').removeClass('hide');
    $('#control-btn-right').removeClass('hide');

    if (dataIdx > 0) {
      $('#control-btn-left').removeClass('disabled');
      $('#control-btn-left:not(disabled)').click(function () {
        my.vis(selector, dataIdx-1);
      });
    } else {
      $('#control-btn-left').addClass('disabled');
    }

    if (dataIdx < my.data.length - 1) {
      $('#control-btn-right').removeClass('disabled');
      $('#control-btn-right:not(disabled)').click(function () {
        my.vis(selector, dataIdx+1);
      });
    } else {
      $('#control-btn-right').addClass('disabled');
    }

    if ($('#checkbox-scale').prop('checked')) {
      my.outerWidth = Math.min($(window).width(), $(window).height()) * 0.8;
    } else {
      my.outerWidth = 600;
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

    var svg = d3.select(selector).append('svg')
      .attr('width', my.outerWidth)
      .attr('height', my.outerHeight);

    var g = svg.append('g')
      .attr('transform', 'translate(' + my.margin.left + ', ' + my.margin.top + ')');

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

    g.selectAll('line.matches')
      .data(data.matches.fwd.concat(data.matches.rev))
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
  };

  return my;
}();

$(maze.main.bind(null, '#vis'));
