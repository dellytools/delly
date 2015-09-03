# suave: structural variant explorer

`suave` is an interactive web application
to visualize read depth ratios
between two samples and the structural variants of one of
the samples (typically the "case" sample in a case/control
setup such as tumor/normal comparison).

## Initial application setup

Firt off, install all required python packages.
You probably want to use `virtualenv` for this.
For example:
```bash
$ virtualenv venv
$ source venv/bin/activate
$ pip install -r requirements/python2.7.txt
```

## Preparing input files

The read depth ratios are computed from HDF5 files, which
in turn are generated from BAM files. In the following example,
we have a BAM file `input.bam` (note: you will need a `.bai`
index file as well) and we generate an output HDF5 file called
`output.h5`. Furthermore, we label this particular sample as
`control` which will make it easier later on to refer to this file. We also apply gzip compression (this is optional but recommended).
```bash
$ python suave_bam_to_h5.py -s control -c gzip -o output.h5 input.bam
```

Furthermore, we need a VCF file. You should use a file that
has been filtered to a set of case-specific calls.

## Starting the server

Now, we can start the server, specifiying one VCF file
and one or more HDF5 files.
By default, the server will listen to port 5000, but
you can change this with the `--port` option.

```bash
$ python suave_server.py -f control.h5 -f tumor01.h5 -f tumor02.h5 -v calls.vcf.gz
```
## Using the GUI

Open a web browser and point it to http://localhost:5000. Click
the setup icon on the top left and specify a case and control sample
together with a chromosome you want to inspect.

The interactive plot is divided into three logical areas:
- the focus track (middle) displays the read depth ratios and let's you zoom and pan across the chromosome
- the context track (bottom) highlights which region is currently in focus. Furthermore, the highlight area can be resized and dragged. By clicking outside of the highlight area, the focus view gets reset to the whole chromosome
- the SV track (top) visualizes the SVs from the VCF file. This track is synchronized with the focus track. SVs that are completely contained in the current viewport are shown as arcs, SVs with only one breakpoint in the visible area are shown as short, vertical lines. You can click on the arcs/lines to snap the viewport to the start and end coordinates of that particular SV, hovering over the arcs/lines displays a tooltip with some information about that SV.
