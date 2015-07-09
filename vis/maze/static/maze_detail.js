// 2014-2015 Markus Hsi-Yang Fritz and Sascha Meiers


// Todo(meiers): Remember to reset network.http.use-cache in Firefox after devlopment

var maze_detail = function () {
  var my = {};

  my.ref   = null;
  my.query = null;
  my.MUMmatches  = null;
  my.LASTmatches = null;
  my.LASTbps     = null;
  my.LASTrefbps  = null;

  my.main = function (visSelector, matchSelector, bpsSelector, refbpsSelector) {

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

    // Get query, ref and MUMer matches from the maze overview page
    if (window['transferData'] != undefined) {
      dat = window.transferData;
      my.ref = dat.ref;
      my.query = dat.query;
      my.MUMmatches = dat.data;
      $('#title-upper-left').append(my.query.name)
    } else {
      $('.spinner').toggleClass('hide', true);
      $(visSelector).append('<div class="alert alert-danger alert-dismissable" role="alert">Received no data from main window.</div>')
      return null;
    }

    // plot dots
    $(visSelector).parent().find(".spinner").toggleClass('hide', true);
    $(visSelector).empty();
    my.vis(visSelector);

    // Get LAST matches and breakpoints
    $.post('/breakpoints', {
          'ref': JSON.stringify(my.ref.seq),
          'query': JSON.stringify(my.query.seq)

        // when successfull
        }, function (result) {

          // matches
          my.LASTmatches = result.matches;
          $(matchSelector).parent().find('.spinner').toggleClass('hide', true);
          my.addLASTmatches(visSelector, matchSelector);
          // breakpoints
          my.LASTbps = result.breakpoints;
          $(bpsSelector).parent().find('.spinner').toggleClass('hide', true);
          my.addLASTbps(bpsSelector);
          // ref breakpoints
          my.LASTrefbps = result.refbreakpoints;
          $(refbpsSelector).parent().find('.spinner').toggleClass('hide', true);
          my.addLASTrefbps(refbpsSelector);
          // event handler
          $('a.list-group-item').mouseover(function(x) {my.toggleMatch(visSelector, $(this).attr('match_index'), 'show')});
          $(' a.list-group-item').mouseout(function(x) {my.toggleMatch(visSelector, $(this).attr('match_index'), 'hide')});

        // when failed
        }, 'json').fail(function(){
          $('.spinner').toggleClass('hide', true);
          $(matchSelector).append('<div class="alert alert-danger alert-dismissable" role="alert">Could not generate last-split matches.</div>');
          $(matchSelector).append('<div class="alert alert-danger alert-dismissable" role="alert">Could not generate breakpoints.</div>');
        });
  };


  // Todo(meiers): oursource duplicate code
  my.vis = function (visSelector) {
    var data = my.MUMmatches;
    var l1 = my.ref.seq.length;
    var l2 = my.query.seq.length;


    $(visSelector).empty();    
    my.outerWidth = $(visSelector).parent().width() * 0.98;
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

    var svg = d3.select(visSelector).append('svg')
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
  my.toggleMatch = function(visSelector, idx, mode) {
    idxs = idx.split(',');
    for (i in idxs) {
      if (mode=='show') $(visSelector + " .geom[match_index=" + idxs[i] + "]").show()
      else              $(visSelector + " .geom[match_index=" + idxs[i] + "]").hide()
    }
  };

  // adds highlightable LAST matches to the dotplot
  my.addLASTmatches = function (visSelector, matchSelector) {
    // draw rectangles
    d3.select(visSelector).select("svg g g.plotarea").selectAll("rect.geom.last")
      .data(my.LASTmatches)
    .enter().append("rect")
      .attr("match_index", function(d,i)   {return i; })
      .attr("x", function(d)      {return my.scales.x(d.d1); })
      .attr("y", function(d)      {return my.scales.y(d.q1); })
      .attr("width", function(d)  {return my.scales.x(d.d2) - my.scales.x(d.d1); })
      .attr("height", function(d) {return my.scales.y(d.q2) - my.scales.y(d.q1); })
      .attr("class", function (d) {return "geom last " + (d.strand == '+' ? 'plus' : 'minus') });
    // list matches (unique id listLASTmatches)
    for (x in my.LASTmatches) {
      var m = my.LASTmatches[x];
      $(matchSelector).append(
            '<a class="list-group-item collapse-group" match_index="' + x + '">' + 
            '  <button class="btn btn-default pull-right" type="button" data-toggle="collapse" data-target="#listLASTmatches a[match_index=' + x + '] pre">' +
            '    <span class="glyphicon glyphicon-collapse-down"></span>' +
            '  </button>' +
            '  <h5>' + 
                 m.sim + ' mismatches. Ref[' + (m.d1) + ':' + (m.d2)  + '] vs. Query[' + m.q1 + ':' + m.q2 + '] (' + m.strand + ')' + 
            '  </h5>' +
            '  <pre class="collapse" aria-expanded="false">' + m.record + '</pre>' +
            '</a>');
    }
    // initially hide rectangles
    $(visSelector + ' svg rect.geom.last').hide();
  };


  my.addLASTbps = function (bpsSelector) {
    for ( x in my.LASTbps) {
      var br = my.LASTbps[x];
      $(bpsSelector).append(
        '<a class="list-group-item collapse-group" match_index="' + br.match_index + '">' + 
            '  <button class="btn btn-default pull-right" type="button" data-toggle="collapse" data-target="' + bpsSelector + ' a[match_index=\'' + br.match_index  + '\'] pre">' +
            '    <span class="glyphicon glyphicon-collapse-down"></span>' +
            '  </button>' +
            '  <h5>' + 
                 (br.length>0 ? br.length + ' bp ' + br.event : 'Exact breakpoint') + ' around query[' + br.qu_coords[0] + ':' + br.qu_coords[1] + ']' +
            '  </h5>' +
            '  <pre class="collapse" aria-expanded="false">' + br.html + '</pre>' +
            '</a>');
    }
  };

  my.addLASTrefbps = function (refbpsSelector) {
    for ( x in my.LASTrefbps) {
      var br = my.LASTrefbps[x];
      console.log(br);
      $(refbpsSelector).append(
        '<a class="list-group-item collapse-group" match_index="' + br.match_index + '">' + 
            '  <button class="btn btn-default pull-right" type="button" data-toggle="collapse" data-target="' + refbpsSelector + ' a[match_index=\'' + br.match_index  + '\'] pre">' +
            '    <span class="glyphicon glyphicon-collapse-down"></span>' +
            '  </button>' +
            '  <h5>' + 
                 (br.length>0 ? br.length + ' bp ' + br.event : 'Exact breakpoint') + ' around Ref[' + br.ref_coords[0] + ':' + br.ref_coords[1] + ']' +
            '  </h5>' +
            '  <pre class="collapse" aria-expanded="false">' + br.html + '</pre>' +
            '</a>');
    }
  };

  return my;
}();

  



/*
for ( x in d.breakpoints) {
                    var br = d.breakpoints[x];
                    $('#div_last_breakpoints ul').append("<li match_index='" + br.match_index + 
                                              "'>Query coordinates: " + br.qu_coords[0] + "/" + 
                                              br.qu_coords[1] + "<br>" +
                                              (br.length>0 ? br.event + ' of length ' + br.length : 'exact break') + 
                                              "<pre>" + br.html + "</pre>");
                  }
*/


// load maze_detail when ready
$(document).ready(function() {
  $(maze_detail.main.bind(null, '#vis', '#listLASTmatches', '#listLASTbps', '#listLASTrefbps'));
});

