DELLY
=====

DELLY is an integrated structural variant prediction method that can detect deletions, tandem duplications, inversions and translocations
at single-nucleotide resolution in short-read massively parallel sequencing data. It uses paired-ends and split-reads to sensitively and accurately
delineate genomic rearrangements throughout the genome.


DELLY dependencies
------------------

* BamTools, (https://github.com/pezmaster31/bamtools)
* Boost C++ Libraries, (www.boost.org)
* zlib compression library, (www.zlib.net)
* kseq library to parse FASTA/FASTQ, (http://lh3lh3.users.sourceforge.net/parsefastq.shtml)

Installing DELLY
----------------

The easiest way to get DELLY is to clone the git repository:

`git clone https://github.com/tobiasrausch/delly.git`

Then you need to install Boost, BamTools and KSEQ and update the paths in the Makefile.
Building DELLY requires just:

`make -B src/delly`

Alternatively, statically linked binaries for Linux 64-bit are available here:

Linux 64-bit Binaries, [delly_v0.0.11.tar.gz](http://www.embl.de/~rausch/delly_v0.0.11.tar.gz)


Running DELLY
-------------


Citation
--------

Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel.

[DELLY: structural variant discovery by integrated paired-end and split-read analysis.](http://bioinformatics.oxfordjournals.org/content/28/18/i333.abstract)

Bioinformatics 2012 28: i333-i339.






<ul>
<li>
Programming Language: C++
</li>
<li>
Supported Platforms: Linux
</li>
<li>
License: Open Source (GNU General Public License, version 3.0)
</li>
</ul>

Please subscribe to the <a href="http://www.embl.de/~rausch/subscribe.html">DELLY mailing list</a> for questions, updates or bug reports.<br>
Previous versions can be downloaded in the <a href="#version">version history</a> below.<br>
<br>


<h3 id="readme" style="background: #F1A340">README</h3>
The DELLY package contains four methods to detect deletions, tandem duplications, inversions and translocations.
There is a separate tool for each SV type:
<ul>
<li  style="list-style-type: square">
Deletions: <b>./delly</b>
</li>
<li  style="list-style-type: square">
Tandem duplications: <b>./duppy</b>
</li>
<li  style="list-style-type: square">
Inversions: <b>./invy</b>
</li>
<li  style="list-style-type: square">
Translocations: <b>./jumpy</b>
</li>
</ul>

Each method has a default help screen with all the command line options. You can simply access it using the method's name, i.e. ./delly, without any additional command line options. Below are typical usage scenarios:
<ul>
<li  style="list-style-type: square">
<b>Paired-end calls: ./delly lib1.bam lib2.bam libN.bam</b><br>
Each bam file is assumed to be a separate library or lane on the sequencer. bam needs to be in lower case otherwise the sam format is assumed.
</li>
<li  style="list-style-type: square">
<b>Paired-end and split-read calls: ./delly -p -g ref.fa file1.bam file2.bam fileN.bam</b><br>
The option '-g' specifies a reference genome. This can be a multi-fasta file. Each header line in the FASTA file has to be a chromosome name identical to the one in the bam file, including all alphanumeric, space or special characters used after the '>' sign in the FASTA file. Many read mappers such as BWA tend to soft-clip one end of the read if it overlaps a structural variant breakpoint. Likewise other read mappers such as ELAND tend to map such reads with sloppy ends containing many mismatches. These reads can be included in the split-read search if you specify the '-p' option, otherwise only single-anchored reads are used.
</li>
<li  style="list-style-type: square">
<b>Paired-end calls using uniquely mapping reads only: ./delly -q 1 -g ref.fa file1.bam file2.bam fileN.bam</b><br>
The '-q' option sets a threshold for the minimum required paired-end mapping quality. All paired-ends with a mapping quality below this threshold are discarded.
</li>
<li  style="list-style-type: square">
<b>Split-read refinement of coarse-grained SV predictions: ./delly -p -e 0.5 -d intervals.txt -g ref.fa file1.bam file2.bam fileN.bam</b><br>
The '-d' option allows to bypass the paired-end step and to directly start the split-read search on the user-provided intervals (intervals.txt). The format of this file is tab-delimited chr, SV start, SV end, Custom SV id. For translocations, it is chr1, pos1, chr2, pos2, id. Use a large value for '-e' (e.g., -e 0.5) if you do expect that your presumed SV intervals are fairly inaccurate.
</li>
<li  style="list-style-type: square">
<b>Parallelizing DELLY / Running DELLY in the cluster: ./delly -p -i your_sample_id -g ref.fa file.bam</b><br>
If you run multiple instances of DELLY in parallel or split your data by chromosome make sure that you prefix all output files by a unique sample ID, e.g.:<br>  
./invy -p -g ref.fa -i ${SAMPLEID} -b ${SAMPLEID}.inv.bp.txt -k ${SAMPLEID}.inv.bp.merge.txt -o ${SAMPLEID}.inv.pe.txt -r ${SAMPLEID}.inv.pe.merge.txt alignment.bam<br>
</li>
<li  style="list-style-type: square">
<b>Excluding certain chromosomes from paired-end and split-read analysis: ./delly -x excludeList.txt -g ref.fa file.bam</b><br>
The '-x' option points to a file with a list of chromosomes that are excluded from paired-end and split-read analysis.
By means of this option you can, for instance, exclude unplaced contigs in your reference genome from structural variant analysis.
The format of the file is a simple list of chromosome names matching the chromosome names in the bam file, e.g.:<br>
chrM<br>
GL000225.1<br>
...
</li>
</ul>

DELLY produces two output files for deletion and tandem duplication calling. One file with the paired-end predictions and one file with the split-read predictions. The split-read predictions are always a subset of the paired-end predictions and the unique structural variant id is kept across the paired-end and split-read calls, so every split-read call can be uniquely linked to a paired-end call.
<ul>
<li  style="list-style-type: square">
<b>Output format of paired-end predictions</b><br>
Each paired-end prediction starts with a print out of all the supporting pairs. The last line is the so-called summary line.
The summary line for the different types of structural variants is tab-delimited and explained below for each SV type, including one example summary line:<br>
<ul>
<li>
<b>Deletion</b><br>
chr, start, end, size, #supporting_pairs, avg._mapping_quality, deletion_id<br>
chr1, 10180, 10509, 329, 75, 15.8667, Deletion_Sample_00000000<br>
</li>
<li>
<b>Tandem Duplication</b><br>
chr, start, end, size, #supporting_pairs, avg._mapping_quality, tandem_duplication_id<br>
chr1, 147772526, 148153764, 381238, 5, 12, Duplication_Sample_00002375<br>
</li>
<li>
<b>Inversion</b><br>
chr, start, end, size, #supporting_pairs, avg._mapping_quality, inversion_id<br>
chr1, 77995222, 106445083, 28449861, 2, 0, Inversion_3to3_Sample_00003084<br>
The inversion id encodes whether the paired-ends are left- (3' to 3') or right-spanning (5' to 5'), please see the DELLY paper for details.<br>
<u>Merged paired-end inversion calls</u><br>
chr, start, end, size, overlap_of_merged_left/right_inversion_calls, avg._mapping_quality_of_both_calls, merged_inversion_ids<br>
chr8, 3769030, 130762835, 126993805, 0.984975, 44, Inversion_5to5_Sample_00061679|Inversion_3to3_Sample_00061732<br>
</li>
<li>
<b>Translocation</b><br>
chrA, posA, chrB, posB, #supporting_pairs, avg._mapping_quality, translocation_id<br>
chr1, 107274529, chr10, 88941781, 4, 65.75, Translocation_3to5_Sample_00003628<br>
The translocation id encodes whether the two chromosomes are joined 5' to 5' (5to5), 3' to 3' (3to3), 3' to 5' (3to5) or 5' to 3' (5to3), please see the DELLY paper for further details.<br>
<u>Merged paired-end translocation calls (translocated segments)</u><br>
chrA, startA, endA, sizeA, chrB/insertion_point/inverted(0=no,1=yes), avg._mapping_quality_of_both_calls/avg._paired-end_support, merged_translocation_ids<br>
chr2, 133026732, 133030695, 3963, chr10/98561434/0, 37/4, Translocation_3to5_Sample_00000121|Translocation_5to3_Sample_00000122
</li>
</ul>
</li>
<li  style="list-style-type: square">
<b>Output format of split-read predictions</b><br>
Each split-read prediction starts with a print out of the consensus alignment of all found split reads.
Afterwards is again a summary line, which follows the format of the paired-end summary line with two exceptions.
First, the number of paired-ends is now replaced by the number of split-reads and second, the average mapping quality is replaced
with the percent alignment quality of the consensus sequence against the reference. Note that the call ids are reused from the paired-end
file so that each split-read call has a corresponding paired-end call. For split-read supported calls the structural variant id is followed by the sequence of the left and right breakpoint, where upper case letters are supported by the read and lower case letters are the next nucleotides in the reference but not in the read (because of the split). The last 4 columns are the length of the microinsertion and the actual microinsertion sequence and the same for any microhomology if present at the breakpoint.
</li>
</ul>
To get a list of all candidate sturctural variants a simple grep for the summary line is sufficient, i.e. grep ">Deletion" del.txt.


<h3 id="faq" style="background: #F1A340">FAQ</h3>
<ul>
<li  style="list-style-type: square">
<b>What is the smallest SV size DELLY can call?</b><br>
This depends on the sharpness of the insert size distribution. For an insert size of 200-300bp with a 20-30bp standard deviation, DELLY starts to call reliable SVs >=150bp.
</li>
<li  style="list-style-type: square">
<b>Can DELLY be used on a non-diploid genome?</b><br>
Yes. At the moment DELLY makes no assumption regarding ploidy of the sample.
</li>
<li  style="list-style-type: square">
<b>What are good cutoffs for the paired-end support and the avg. paired-end mapping quality?</b><br>
This depends on many factors, including the sequencing coverage, the repetitiveness of the genome (and the loci of interest), the insert size distribution, and so on. Taking these caveats into account a good start might be a paired-end cutoff >=3 and a mapping quality >=20.
</li>
<li  style="list-style-type: square">
<b>Is there any visualization of the called SVs?</b><br>
No, DELLY does not produce any graphical output. However, there are many viewers such as the <a href="http://www.broadinstitute.org/igv/">IGV</a> that do visualize discordantly mapped paired-ends. Thus, localizing and investigating these events in the IGV browser indeed makes sense.
</li>
<li  style="list-style-type: square">
<b>How can DELLY be used to call somatic SVs?</b><br>
Run DELLY separately on as many cancer and control genomes you have access to and then filter the tumor SVs of a single sample against all SVs called in all control genomes. In practice, this appears to be one of the best options to derive a high quality set of somatic SVs. For copy-number variable events (CNVs) such as deletions and tandem duplications a further annotation with read-depth may be helpful, although most complex rearrangements do not necessarily show such a read-depth change.
</li>
<li  style="list-style-type: square">
<b>What directory is used by DELLY to store temporary files?</b><br>
DELLY uses the directory in the TMP environment variable for temporary files.
</li>
<li  style="list-style-type: square">
<b>Can DELLY be used on bwa mem alignments?</b><br>
DELLY can be used with bwa mem but you have to mark shorter split-read alignments as secondary alignments using the '-M' option, e.g. bwa mem -M ref.fa file.fq.gz.
</li>
<li  style="list-style-type: square">
<b>Are non-unique alignments, multi-mappings and/or multiple split-read alignments allowed?</b><br>
DELLY expects two alignment records in the bam file for every paired-end, one for the first and one for the second read. Multiple split-read alignment records of a given read are allowed if and only if one of them (e.g. the longest split alignment) is a primary alignment whereas all others are marked as secondary (flag 0x0100).
</li>
</ul>

<h3 id="version" style="background: #F1A340">Version History</h3>
<ul>
<li  style="list-style-type: square">
<b>v0.0.11:</b><a href="delly_v0.0.11.tar.gz">Binaries</a>, <a href="delly_source_v0.0.11.tar.gz">Source Code</a> (6th June 2013)<br>
<ul>
<li>bwa mem is now supported if and only if the bwa mem '-M' option is used to mark short split-read alignments as secondary alignments leaving only one primary alignment for any given read.</li>
<li>A split-read alignment progress bar has been added.</li>
<li>Inversions are now annotated by their configuration using 3' to 3' (3to3) for lef-spanning inversions and 5' to 5' for right-spanning inversions.</li>
</ul></li>
<li  style="list-style-type: square">
<b>v0.0.10:</b><a href="delly_v0.0.10.tar.gz">Binaries</a>, <a href="delly_source_v0.0.10.tar.gz">Source Code</a> (23rd May 2013)<br>
<ul>
<li>Fasta identifiers are truncated at the first whitespace.</li>
<li>Unique temporary filenames, no need to specify '-i' anymore.</li>
<li>Dropped '-z'. All temporary files are written to the directory specified in the TMP environment variable.</li>
<li>Cluster paired-end mapping quality is now the median instead of the mean.</li>
</ul></li>
<li  style="list-style-type: square">
<b>v0.0.9:</b><a href="delly_v0.0.9.tar.gz">Binaries</a>, <a href="delly_source_v0.0.9.tar.gz">Source Code</a> (1st November 2012)<br>
<ul>
<li>Unified the split-read search across all tools.</li> 
<li>Dropped 1MB size restriction for split-read refinement of inversions, deletions and tandem duplications.</li>
<li>Reference index uses now k=11 instead of k=7 to slightly increase specificity in the split-read search (and because reads got longer >75bp).</li>
<li>For invy, paired-end SV intervals are now cut in the middle and then one side is flipped.</li>
<li>The average quality computation now assumes Phred quality scores.</li>
<li>Minor changes in the output format of DELLY.</li>
</ul></li>
<li  style="list-style-type: square">
<b>v0.0.8:</b><a href="delly_v0.0.8.tar.gz">Binaries</a>, <a href="delly_source_v0.0.8.tar.gz">Source Code</a> (12th October 2012)<br>
<ul>
<li>New '-x' option to exclude a list of chromosomes from PE and SR analysis.</li>
<li>Replaced vertex->sam record hash map with direct pointers.</li>
</ul></li>
<li  style="list-style-type: square">
<b>v0.0.7:</b><a href="delly_v0.0.7.tar.gz">Binaries</a>, <a href="delly_source_v0.0.7.tar.gz">Source Code</a> (7th September 2012)<br>
<ul>
<li>Split-read positions are 1-based.</li>
<li>Split-read SV start position is the last base before the SV.</li>
<li>Split-read SV end position is the last base of the SV. Hence, end - start is the SV size.</li>
<li>Reduced RAM usage since nodes are now identified by pointers and not read names anymore.</li>
</ul></li>
<li  style="list-style-type: square">
<b>v0.0.6:</b><a href="delly_v0.0.6.tar.gz">Binaries</a>, <a href="delly_source_v0.0.6.tar.gz">Source Code</a> (3rd September 2012)<br>
<ul>
<li>Libraries/Input bam files are now tagged as Library0, Library1, etc.</li>
<li>In the paired-end output, supporting paired-ends are labeled by this new library tag.</li>
<li>All tools can now be run in split-read mode only, a small example was added to the readme above.</li>
<li>The inverted/non-inverted column was removed from jumpy's output files.</li>
<li>Translocations are further annotated now by their configuration, either 5' to 5' (5to5), 3' to 5' (3to5), 5' to 3' (5to3) or 3' to 3' (3to3).</li>
</ul></li>
<li  style="list-style-type: square">
<b>v0.0.5:</b><a href="delly_v0.0.5.tar.gz">Binaries</a>, <a href="delly_source_v0.0.5.tar.gz">Source Code</a> (26th July 2012)<br>
<ul>
<li>Fixed a minor bug in jumpy's split-read search for translocations starting at zero.</li>
<li>Translocation type (inOrder - yes/no, inverted - yes/no) is now encoded in the translocation id.</li>
<li>Added final date/time stamp at the end of the program.</li>
</ul></li>
<li  style="list-style-type: square">
<b>v0.0.4:</b><a href="delly_v0.0.4.tar.gz">Binaries</a>, <a href="delly_source_v0.0.4.tar.gz">Source Code</a> (24th July 2012)<br>
<ul>
<li>SD calculation was replaced by the Median Absolute Deviation (MAD).</li>
<li>For deletion calling, new parameter '-c' to specify directly a fixed fraction of discordant PEs used for variant calling. This option is experimental, the default MAD is recommended.</li>
<li>Added date/time stamp and command-line parameters to stdout for logging.</li>
</ul></li>
<li  style="list-style-type: square">
<b>v0.0.3:</b><a href="delly_v0.0.3.tar.gz">Binaries</a>, <a href="delly_source_v0.0.3.tar.gz">Source Code</a> (1st July 2012)<br>
<ul>
<li>Code was updated to support the latest releases of BamTools and Boost.</li>
<li>PE clustering graph uses less memory now.</li>
<li>Dropped two unnecessary maps for better performance.</li>
<li>Modified the command line options:</li>
<li>*Dropped '-c' and '-q', always all PE calls are passed to SR refinement now.</li>
<li>*Reused '-q' to specify a min. paired-end mapping quality cutoff.</li>
<li>*Added '-z' to the generic options for the tmp path.</li>
<li>*Dropped '-n' and '-f'. More stringent cutoffs for the #split-reads or the alignment quality can be simply applied to the output of DELLY (SR).</li>
</ul></li>
<li  style="list-style-type: square">
<b>v0.0.2:</b><a href="delly_v0.0.2.tar.gz">Binaries</a>, <a href="delly_source_v0.0.2.tar.gz">Source Code</a> (11th June 2012)<br>
<ul>
<li>Lanes from the same library can now be merged in a single large bam file.</li>
<li>Optimized removal of redundant pairs.</li>
</ul></li>
<li  style="list-style-type: square">
<b>v0.0.1:</b><a href="delly_v0.0.1.tar.gz">Binaries</a>, <a href="delly_source_v0.0.1.tar.gz">Source Code</a> (22nd March 2012)<br>
<ul>
<li>Initial release of DELLY.</li>
</ul></li>
</ul>


<h3 id="links" style="background: #F1A340">Links</h3>
<ul>
<li  style="list-style-type: square">
<a href="http://www.embl.de/">EMBL</a>
</li>
<li  style="list-style-type: square">
<a href="http://www.korbel.embl.de">Korbel Labpage</a>
</li>
<li  style="list-style-type: square">
<a href="http://genecore3.genecore.embl.de/genecore3/">EMBL GeneCore</a>
</li>
<li  style="list-style-type: square">
<a href="http://www.embl.de/~rausch/">Personal Homepage</a>
</li>
</ul>

<h3 id="acknowledgements" style="background: #F1A340">Acknowledgements</h3>
<ul>
<li  style="list-style-type: square">
Boost C++ Libraries, <a href="http://www.boost.org/">Boost</a>
</li>
<li  style="list-style-type: square">
BamTools (Derek Barnett), <a href="http://sourceforge.net/projects/bamtools/">BamTools</a>
</li>
<li  style="list-style-type: square">
zlib compression library, <a href="http://zlib.net">zlib</a>
</li>
<li  style="list-style-type: square">
SeqAn library, <a href="http://www.seqan.de/">SeqAn</a>
</li>
<li  style="list-style-type: square">
1000 Genomes Project, <a href="http://www.1000genomes.org/">1000GP</a>
</li>
<li  style="list-style-type: square">
International Cancer Genome Consortium, <a href="http://www.icgc.org/">ICGC</a>
</li>
</ul>
<hr>
<hr>


<table border="0" cellspacing="20" width="100%">
  <tr>
    <td>
    <p>
    <a href="http://validator.w3.org/check?uri=referer"><img border=0 src="http://www.w3.org/Icons/valid-html401" alt="Valid HTML 4.01 Transitional" height="31" width="88"></a>
    </p>
    <p>
    <a href="http://jigsaw.w3.org/css-validator/check/referer"><img style="border:0;width:88px;height:31px" src="http://jigsaw.w3.org/css-validator/images/vcss" alt="CSS ist valide!"></a>
    </p>
    </td>
    <td>
    <p align="right">
    <small><a href="http://www.embl.de/~rausch/disclaimer.html">*EMBL Disclaimer</a></small>
    <br>
    Author: <a href="http://www.embl.de/~rausch/">Tobias Rausch</a>
    <br>
    Last change: 6th June 2013
    <br>
    Main Page of the <a href="http://www.embl-heidelberg.de/" target="_top">EMBL-Heidelberg</a>
    </p>
    </td>
  </tr>
</table>

</body></html>

