// 2014-2015 Markus Hsi-Yang Fritz and Sascha Meiers


// Todo(meiers): Remember to reset network.http.use-cache in Firefox after devlopment

var maze_detail = function () {
  var my = {};

  my.ref   = null;
  my.query = null;
  my.MUMmatches  = null;
  my.LASTmatches = null;

  my.main = function (selector) {

    // Todo(meiers): outsource duplicate code
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

    // Get query, ref and matches from the maze overview page
    if (window['transferData'] != undefined) {
      dat = window.transferData;
      my.ref = dat.ref;
      my.query = dat.query;
      my.MUMmatches = dat.data;
    } else {
      alert('FATAL: Did not receive data from main window.')
      return null;
    }

    // plot dots
    $('#title-upper-left').append(my.query.name)
    $(selector).empty();
    $('.spinner').toggleClass('hide');
    my.vis(selector);

    // get LAST matches
    $.post('/breakpoints', {
          'ref': JSON.stringify(my.ref.seq),
          'query': JSON.stringify(my.query.seq)
        }, function (res) {
          $('.spinner2').toggleClass('hide');
          my.LASTmatches = res.matches;
          
          my.addLASTmatches(selector);
        }, 'json'
      );
  };


  // Todo(meiers): oursource duplicate code
  my.vis = function (selector) {
    var data = my.MUMmatches;
    var l1 = my.ref.seq.length;
    var l2 = my.query.seq.length;

    $(selector).empty();    
    my.outerWidth = Math.min(500, $(window).height()) * 0.8;
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
      .attr('clip-path', 'url(#clipChartArea)')
      .attr('class', 'plotarea');

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
        /*
        chartArea.selectAll('.geom')
          .data(data.fwd.concat(data.rev))
          .attr('x1', function (d) { return x(d[0]); })
          .attr('x2', function (d) { return x(d[1]); })
          .attr('y1', function (d) { return y(d[2]); })
          .attr('y2', function (d) { return y(d[3]); });
        */
        chartArea.selectAll(".geom")
           .attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
        g.select(".x.axis.b").call(xAxisB);
        g.select(".x.axis.t").call(xAxisT);
        g.select(".y.axis.r").call(yAxisR);
        g.select(".y.axis.l").call(yAxisL);
      }

      // make scales accessible from outside:
      my.scales = { x: x, y: y}
  };

  // show or hide all elements belonging to a match
  my.toggleMatch = function(selector, idx, mode='show') {
    console.log("Toggle: " + idx)
    idxs = idx.split(',');
    for (i in idxs) {
      console.log(idxs[i])
      if (mode=='show') $(selector + " .geom[match_index=" + idxs[i] + "]").show()
      else              $(selector + " .geom[match_index=" + idxs[i] + "]").hide()
    }
  };

  // Todo(meiers): consider use of selector here. Better use 2 selectors 
  //              (one for plot, one for output table).
  // adds highlightable LAST matches to the dotplot
  my.addLASTmatches = function (selector) {
    
    // draw rectangles
    d3.select(selector).select("svg g g.plotarea").selectAll("rect.geom.last")
      .data(my.LASTmatches)
    .enter().append("rect")
      .style("opacity", "0.3")
      .attr("class", "geom last")
      .attr("match_index", function(d,i)   {return i; })
      .attr("x", function(d)      {return my.scales.x(d.d1); })
      .attr("y", function(d)      {return my.scales.y(d.q1); })
      .attr("width", function(d)  {return my.scales.x(d.d2) - my.scales.x(d.d1); })
      .attr("height", function(d) {return my.scales.y(d.q2) - my.scales.y(d.q1); })
      .style("fill", function (d) {return d.strand == '+' ? 'black' : 'red'});
      

      // Todo(meiers): CSS define matchlist and matchlistHeader
    // list matches
    $('div#match-wrap').append('<span class="matchlistHeader">last-split matches</span>' + 
                                '<ul class="matchlist" id="listLASTmatches"></ul>');
    for (x in my.LASTmatches) {
      var m = my.LASTmatches[x];
      $('div#match-wrap ul#listLASTmatches').append("<li class='" + (m.strand == '+' ? 'plus' : 'minus') + 
                    "' match_index='" + x + "'>" + m.sim + " mismatches. " + 
                    "Ref[" + (m.d1) + ":" + (m.d2)  + "] (" + m.d1 + ":" + 
                    m.d2 + ") query[" + m.q1 + ":" + m.q2 + "] <pre>" + m.record + 
                    "</pre></li>");
    }

    // Todo(meiers): This event handler will probably hold even for breakpoints etc.
    // Event handler
    $('div#match-wrap li').click(function(d) { $('pre', this).toggle(); });
    $('div#match-wrap li').mouseover(function(x) {
                      my.toggleMatch(selector, $(this).attr('match_index'), 'show')});
    $('div#match-wrap li').mouseout(function(x) {
                      my.toggleMatch(selector, $(this).attr('match_index'), 'hide')});

    // initially hide
    $('#chart svg svg rect.geom.last').hide();
    $('div#match-wrap li pre').hide();
  };

  return my;
}();


// load maze_detail when ready
$(document).ready(function() {
  $(maze_detail.main.bind(null, '#vis'));
});

