#! /usr/bin/env python

from __future__ import print_function
import vcf
import argparse
import numpy
import collections
import gc
import math
import pylab
import matplotlib
import matplotlib.pyplot
import matplotlib.gridspec
import sklearn.mixture
import sklearn.metrics
import sklearn.metrics.pairwise

def cnOffset(rcSV, rcLeft, rcRight, geno):
    hetRC = list()
    refRC = list()
    for sp in geno.keys():
        rcLR = rcLeft[sp] + rcRight[sp]
        if rcLR > 0:
            if sum(geno[sp]) == 0:
                refRC.append(2.0 * float(rcSV[sp]) / float(rcLR))
            elif sum(geno[sp]) == 1:
                hetRC.append(2.0 * float(rcSV[sp]) / float(rcLR))
    if len(hetRC) and len(refRC):
        return abs(numpy.median(numpy.array(hetRC)) - numpy.median(numpy.array(refRC)))
    else:
        return None

# Find GMM
def fitData(counts, inputMeans):
    gmm = sklearn.mixture.GMM(params='wc', init_params='wc', n_iter=0, n_components=len(inputMeans)).fit(inputMeans)
    gmm.set_params(n_iter=100)
    try:
        model = gmm.fit(counts)
    except RuntimeError:
        return (None, 0.0)
    classPredict = model.predict(counts)
    sc = 0.0
    if numpy.unique(classPredict).size > 1:
        points = [[0, x] for x in counts]
        pairwiseDist = sklearn.metrics.pairwise.pairwise_distances(points, points, metric='euclidean')
        sc = sklearn.metrics.silhouette_score(pairwiseDist, classPredict, metric='euclidean')
    return (model, sc)

# Bin resolution
def numberOfBins(data):
    q1 = numpy.percentile(data, 25)
    q3 = numpy.percentile(data, 75)
    n = data.size**(.1/.3)
    rng = numpy.amax(data) - numpy.amin(data)
    iqr = 2*(q3-q1)
    bins = ((n*rng)/iqr)
    if (math.isnan(bins)) or (math.isinf(bins)) or (bins < 50):
        return 75
    else:
        return int(bins)

# Plotting
def plotRD(readCount, comp, model, pTitle, hCount):
    fig = matplotlib.pyplot.figure(figsize=(12, 6))
    gs = matplotlib.gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[2, 1])
    ax = fig.add_subplot(gs[0, :])
    bins = numpy.linspace(min(readCount), max(readCount), numberOfBins(readCount))
    lp, rp = model.score_samples(bins)
    pdf = numpy.exp(lp)
    pdf_individual = rp * pdf[:, numpy.newaxis]
    ax.hist(readCount, bins, normed=True, histtype='stepfilled', alpha=0.4)
    ax.plot(bins, pdf_individual, '-k')
    for clustIndex in range(comp):
        ax.axvline(int(model.means_[clustIndex]), linewidth=1.0, color="darkgreen")
    if hCount is not None:
        ax.axvline(hCount, linewidth=1.0, color="darkred")
    ax.set_xlim(min(bins), max(bins))
    ax.set_xlabel("Normalized read counts", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_title(pTitle, fontsize=10)

    az = fig.add_subplot(gs[1, :])
    p = model.predict_proba(bins)
    cm = pylab.get_cmap('Dark2')
    for index in range(comp):
        colVal = cm(1.*index/comp)
        az.fill_between(bins, p[:, index], color=colVal, alpha=0.5)
    for index in range(comp):
        az.plot(bins, p[:, index], color='black')
    az.set_xlim(min(bins), max(bins))
    az.set_ylim(0, 1)
    az.set_xlabel('Normalized read counts', fontsize=10)
    az.set_ylabel(r'$p({\rm Class}|ReadCount)$', fontsize=10)

    fig.tight_layout()
    fileName = "" + svID + ".png"
    fig.savefig(fileName, dpi=100)
    fig.clf()
    pylab.close()
    matplotlib.pyplot.close('all')



# Parse command line
parser = argparse.ArgumentParser(description='Deletion/Duplication filter.')
parser.add_argument('-v', '--vcf', metavar='cnv.vcf', required=True, dest='cnvVCF', help='deletion/duplication vcf file (required)')
parser.add_argument('-o', '--outVCF', metavar='out.vcf', required=True, dest='outVCF', help='output vcf file (required)')
parser.add_argument('-i', '--id', metavar='esv2659357', required=False, dest='cnvID', help='only compute SV id (optional)')
parser.add_argument('-l', '--label', metavar='NA12878', required=False, dest='highlightID', help='label given sample (optional)')
parser.add_argument('-p', '--plot', dest='doPlot', action='store_true', help='plot read-depth')
args = parser.parse_args()

# Parse id
cnvID = None
if args.cnvID:
    cnvID = str(args.cnvID)

# Parse id
highlightID = None
if args.highlightID:
    highlightID = str(args.highlightID)


# Parse deletions/duplications
if args.cnvVCF:
    vcf_reader = vcf.Reader(open(args.cnvVCF), 'r', compressed=True) if args.cnvVCF.endswith('.gz') else vcf.Reader(open(args.cnvVCF), 'r', compressed=False)
    samples = vcf_reader.samples
    vcf_reader.infos['SCORE'] = vcf.parser._Info('SCORE', 1, 'Float', 'Silhouette score')
    vcf_reader.infos['CNOFFSET'] = vcf.parser._Info('CNOFFSET', 1, 'Float', 'CN Offset')
    vcf_reader.formats['CNL'] = vcf.parser._Format('CNL', '.', 'Float', 'Copy-number genotype likelihoods')
    vcf_writer = vcf.Writer(open(args.outVCF, 'w'), vcf_reader, lineterminator='\n')
    vcf_writer.close()
    with open(args.outVCF, 'a') as f:
        seen = False
        for record in vcf_reader:
            if seen:
                break
            if cnvID is not None:
                if cnvID != record.ID:
                    continue
                else:
                    seen = True
            gt = dict()
            rc = collections.defaultdict(int)
            rcl = collections.defaultdict(int)
            rcr = collections.defaultdict(int)
            peCount = 0
            for call in record.samples:
                rc[call.sample] = call['RC']
                rcl[call.sample] = call['RCL']
                rcr[call.sample] = call['RCR']
                if call.called:
                    gt[call.sample] = [int(gVal) for gVal in call['GT'].split('/')]
                    if (call.gt_type != 0) and (call['DV'] > 0):
                        peCount += call['DV']
            if len(gt):
                cno = cnOffset(rc, rcl, rcr, gt)
                if (cno > 0.5) and (cno < 1.5):
                    chrName = record.CHROM
                    start = record.POS
                    end = record.INFO['END']
                    svID = record.ID

                    # Load read counts
                    count = numpy.arange(len(samples), dtype=numpy.float)
                    ind = 0
                    highlightCount = None
                    for ind, s in enumerate(samples):
                        den = rcl[s] + rcr[s]
                        if den > 0:
                            count[ind] = 2 * float(rc[s])/ (float(den))
                        else:
                            count[ind] = None
                        if (highlightID is not None) and (s == highlightID):
                            highlightCount = count[ind]

                    # Fit the GMM and get the deviation of the fit
                    cnstates = range(9)
                    M_best, score = fitData(count, cnstates)
                    if score == 0.0:
                        continue
                    cncluster = M_best.predict(cnstates)
                    clusterToCN = dict()
                    for ind, cl in enumerate(cncluster):
                        clusterToCN[cl] = ind

                    # Get predicted CN for every sample
                    logprob, responsibilities = M_best.score_samples(count)
                    cluster = M_best.predict(count)
                    cnSample = dict()
                    for ind, s in enumerate(samples):
                        try:
                            cnSample[s] = clusterToCN[cluster[ind]]
                        except KeyError:
                            cnSample[s] = -1

                    # Compute GLs
                    sampleGL0 = dict()
                    sampleFT = dict()
                    for sample in samples:
                        sampleFT[sample] = "LowQual"
                        for i in cnstates:
                            sampleGL0[(sample, i)] = 1e-300
                    for ind, sample in enumerate(samples):
                        for cl in numpy.unique(cluster):
                            if (responsibilities[ind, cl] > 1e-300) and (cl in clusterToCN.keys()):
                                sampleGL0[(sample, clusterToCN[cl])] = responsibilities[ind, cl]
                            if responsibilities[ind, cl] > 0.95:
                                sampleFT[sample] = "PASS"

                    # Output
                    score = numpy.round(score, decimals=2)
                    cno = numpy.round(cno, decimals=2)

                    info = "IMPRECISE;"
                    info += "SVTYPE=CNV;"
                    info += "PE=" + str(record.INFO['PE']) + ";"
                    info += "MAPQ=" + str(record.INFO['MAPQ']) + ";"
                    info += "CHR2=" + str(record.INFO['CHR2']) + ";"
                    info += "END=" + str(record.INFO['END']) + ";"
                    info += "SCORE=" + str(score) + ";"
                    info += "CNOFFSET=" + str(cno) + ";"
                    gtStr = "CNL:FT:RCL:RC:RCR:CN:DR:DV:RR:RV"
                    print(record.CHROM, record.POS, record.ID, record.REF, record.ALT[0], ".", "PASS", info, gtStr, sep="\t", file=f, end="")
                    for call in record.samples:
                        cnl0Str = str(numpy.round(math.log10(sampleGL0[(call.sample, 0)])))
                        for i in cnstates[1:]:
                            cnl0Str += "," + str(numpy.round(math.log10(sampleGL0[(call.sample, i)])))
                        print("\t", file=f, end="")
                        print(cnl0Str, sampleFT[call.sample], call['RCL'], call['RC'], call['RCR'], cnSample[call.sample], call['DR'], call['DV'], call['RR'], call['RV'], sep=':', file=f, end='')
                    print("", file=f)

                    # Plot the clustering
                    if args.doPlot:
                        plotTitle = "" + str(chrName) + ":" + str(start) + "-" + str(end) + "; " + svID + " (SV-Size=" + str(end-start) + ", Silhouette=" + str(score) + ", CN-Offset=" + str(cno) + ")"
                        plotRD(count, len(cnstates), M_best, plotTitle, highlightCount)
                        gc.collect()
