import argparse, os, sys, logging
import random
import shutil
import pysam

parser = argparse.ArgumentParser(description='')

parser.add_argument('all_novel', help='A file containing information for all the novel sites')
parser.add_argument('-a', '--annotation', default='/projects/dmacmillanprj2/polya/ccle/novel_cleavage_events/pysam_version/with_utr3/ucsc.utr3.gtf', help='UCSC-downloaded annotation file')
parser.add_argument("-l", "--log", dest="logLevel", default='WARNING', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help='Set the logging level. Default = "WARNING"')
parser.add_argument('-n', '--name', default='result', help='Name for the file output. Default is "result"')
parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is {}'.format(os.getcwd()))

args = parser.parse_args()

filtered_kleats_path = '/projects/dmacmillanprj2/polya/ccle/filteredKleats/kleats_plus_added'
webdir = '/gsc/www/bcgsc.ca/downloads/dmacmillan/'
data_output = open(os.path.join(args.outdir, 'data.compareNovel'), 'w')

def getGtexCram(name):
    return '/projects/btl/polya/gtex/star/crams/{}.star.cram'.format(name)

def getCcleCram(name):
    return '/projects/btl/polya/sorted_star_ccle/cram_sorted_{}.cram'.format(name)

if not os.path.isdir(args.outdir):
    try:
        os.makedirs(args.outdir)
    except OSError:
        pass

# Logging
logging.basicConfig(filename=os.path.join(args.outdir, 'log.compareNovel'), level=getattr(logging, args.logLevel))

novel = strong = gtex = utrs = None

with open(args.annotation, 'r') as f:
    utrs = [x.strip().split('\t') for x in f.readlines()]

with open(args.all_novel) as f:
    all_novel = [x.strip().split('\t') for x in f.readlines()]
    all_novel = all_novel[1:]

all_ccle_samples = set([x[0] for x in all_novel])

gnovel = {}

for i in all_novel:
    if i[4] not in gnovel:
        gnovel[i[4]] = {}
    if i[5] not in gnovel[i[4]]:
        gnovel[i[4]][i[5]] = []
    gnovel[i[4]][i[5]].append(i)

def sprint(text):
    sys.stdout.write(text)
    sys.stdout.flush()

def getClosestUtr3(site, utrs):
    utrs = [x for x in utrs if x[0] == site[4]]
    min_dist = float('inf')
    closest_utr3 = None
    for i, u in enumerate(utrs):
        start = int(u[3])
        end = int(u[4])
        cs = int(site[6])
        dist = min(abs(start-cs),abs(end-cs))
        if dist < min_dist:
            min_dist = dist
            closest_utr3 = u
    return closest_utr3


window = 20
clusters = []
for chrom in gnovel:
    for gene in gnovel[chrom]:
        cluster = {'centroid': None, 'sites': []}
        for nov in gnovel[chrom][gene]:
            cs = int(nov[6])
            if not cluster['centroid']:
                cluster['centroid'] = cs
                cluster['sites'].append(nov)
                continue
            dist = abs(cluster['centroid'] - cs)
            if dist <= 20:
                cluster['sites'].append(nov)
                cluster['centroid'] = sum([int(x[6]) for x in cluster['sites']]) / len(cluster['sites'])
                continue
            else:
                clusters.append(cluster)
                break

clusters = sorted(clusters, key=lambda x: len(x['sites']), reverse=True)

def getMedian(array):
    length = len(array)
    if (length == 0):
        return 'N/A'
    elif (length % 2 == 0):
        return ((array[length/2]) + (array[(length/2)-1]))/float(2)
    else:
        return array[length/2]

def getMedianCoverage(pysam_alignment_file_object, chrom, start, end):
    ns = []
    for p in pysam_alignment_file_object.pileup(chrom, start, end, truncate=True):
        ns.append(p.n)
    median = getMedian(ns)
    return median

centroids = {}
data_output.write(('\t').join(['CHROM', 'GENE', 'STRAND', 'TISSUE', 'DIST_FROM_ANNOT',
                               'PAS', 'ID', 'CELL_LINE', 'CLEAVAGE_SITE_CENTROID',
                               'SCORE', 'MEDIAN_LEFT', 'MEDIAN_RIGHT',
                               'MED_DIFF', 'CLOSEST_UTR3_ATTR']) + '\n')
for clust in clusters:
    clust['best_site'] = None
    clust['best_score'] = 0
    logging.debug('Total sites in cluster: {}'.format(len(clust['sites'])))
    for i,site in enumerate(clust['sites']):
        logging.debug('site index: {}'.format(i))
        logging.debug('site: {}'.format(site))
        try:
            score = (int(site[16]) + int(site[17]) + int(site[18]) + int(site[19])) / int(site[15][1])
        except (ValueError, IndexError) as e:
            score = (int(site[16]) + int(site[17]) + int(site[18]) + int(site[19]))
        clust['sites'][i].append(score)
        if score > clust['best_score']:
            clust['best_site'] = site
            clust['best_score'] = score
        clust['best_score'] = score
        clust['best_site'] = site
        if not clust['best_site']:
            logging.debug('Skipped because no best site')
            continue
        logging.debug('Best site: {}'.format(clust['best_site']))
        closest_utr3 = getClosestUtr3(clust['best_site'], utrs)
        if not closest_utr3:
            logging.debug('Skipped because no close utr3')
            continue
        logging.debug('Closest utr3: {}'.format(closest_utr3))
        name = clust['best_site'][0]
        cline = clust['best_site'][1]
        annot_dist = clust['best_site'][11]
        pas = clust['best_site'][15]
        tissue = clust['best_site'][2]
        chrom = clust['best_site'][4]
        gene = clust['best_site'][5]
        cs = int(clust['best_site'][6])
        strand = clust['best_site'][10]
        utr3_start = int(closest_utr3[3])
        utr3_end = int(closest_utr3[4])
        if (utr3_start >= cs) or (utr3_end <= cs):
            logging.debug('Skipped because cs not within utr3')
            continue
        #print '-'*(cs-utr3_start) + '|' + '-'*(utr3_end-cs)
        ccle_aln = pysam.AlignmentFile(getCcleCram(name), 'rc')
        left = getMedianCoverage(ccle_aln, chrom, utr3_start, cs)
        right = getMedianCoverage(ccle_aln, chrom, cs, utr3_end)
        if strand == '-':
            temp = right
            right = left
            left = temp
        try:
            diff = left - right
        except TypeError:
            logging.error('TypeError abs(left - right)\nleft: "{}"\nright: "{}"'.format(left, right))
            diff = 0
        logging.debug('diff: {}'.format(diff))
        data = [chrom, gene, strand, tissue, annot_dist, pas, name, cline, clust['centroid'], clust['best_score'], left, right, diff, closest_utr3[8]]
        data_output.write('\t'.join([str(x) for x in data]) + '\n')
