// 2014-2015 Markus Hsi-Yang Fritz

var maze = function () {
  var my = {};

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

    $('#visualize').click(function () {
      $('#configModal').modal('hide');
      $('.spinner').toggleClass('hide');

      var matches = $('#config-matches input:radio:checked').attr('id');
      var length = $('#config-length').val();

      $.getJSON('/data',
        {'matches': matches, 'length': length},
        function (data) {
          $('.spinner').toggleClass('hide');
          console.log(data);
        }
      );
    });
  };

  return my;
}();

$(maze.main.bind(null, '#vis'));
