// Copyright (c) 2014 Markus Hsi-Yang Fritz

var suave = {

  data: null,

  main: function (selector) {
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
      $.getJSON('/data/' + s1 + '/' + s2 + '/' + c, function (res) {
        suave.data = res;
        suave.vis(selector);
      });
    });
  },

  vis: function(selector) {
    $(selector).html('vis()');
  }

};

$(suave.main.bind(null, '#vis'));
