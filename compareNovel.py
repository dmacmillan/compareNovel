import argparse, os, sys, logging
import random
import shutil
import pysam
import urllib
import requests
from bs4 import BeautifulSoup

parser = argparse.ArgumentParser(description='')

parser.add_argument('all_novel', help='A file containing information for all the novel sites')
parser.add_argument('strong_novel', help='A file containing novel sites and their relative frequencies')
parser.add_argument('gtex_list', help='File containing information on selected GTEX samples')
parser.add_argument('kleats_filtered', nargs='+', help='Filtered KLEAT files')
parser.add_argument('-t', '--threshold', type=int, default=100, help='The minimum difference between medians of 3\'UTR regions (as divided by novel site) to be reported. Default = 100')
parser.add_argument('-a', '--annotation', default='/projects/dmacmillanprj2/polya/ccle/novel_cleavage_events/pysam_version/with_utr3/ucsc.utr3.gtf', help='UCSC-downloaded annotation file')
parser.add_argument("-l", "--log", dest="logLevel", default='WARNING', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help='Set the logging level. Default = "WARNING"')
parser.add_argument('-n', '--name', default='result', help='Name for the file output. Default is "result"')
parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is {}'.format(os.getcwd()))

args = parser.parse_args()

filtered_kleats_path = '/projects/dmacmillanprj2/polya/ccle/filteredKleats/kleats_plus_added'
webdir = '/gsc/www/bcgsc.ca/downloads/dmacmillan/'

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
logging.basicConfig(filename=os.path.join(args.outdir, 'compareNovel.log'), level=getattr(logging, args.logLevel))

novel = strong = gtex = utrs = None

with open(args.annotation, 'r') as f:
    utrs = [x.strip().split('\t') for x in f.readlines()]

with open(args.gtex_list, 'r') as f:
    all_gtex = [x.strip().split('\t') for x in f.readlines()]

with open(args.all_novel) as f:
    all_novel = [x.strip().split('\t') for x in f.readlines()]
    all_novel = all_novel[1:]

gnovel = {}

for i in all_novel:
    if i[4] not in gnovel:
        gnovel[i[4]] = {}
    if i[5] not in gnovel[i[4]]:
        gnovel[i[4]][i[5]] = []
    gnovel[i[4]][i[5]].append(i)

with open(args.strong_novel) as f:
    strong = [x.strip().split('\t') for x in f.readlines()]
    #strong = strong[1:]

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

def detectSwitch(site, all_gtex, utrs, gnovel):
    logging.debug('Detecting switch for {} ...'.format(site))
    best_site = compare(site, gnovel)
    logging.debug('Best site: {}'.format(best_site))
    chrom = best_site[4]
    closest_utr3 = getClosestUtr3(best_site, utrs)
    logging.debug('Closest UTR3: {}'.format(closest_utr3))
    #print closest_utr3
    cs = int(best_site[6])
    logging.debug('Novel site from best_match: {}'.format(cs))
    if not (int(closest_utr3[3]) < cs < int(closest_utr3[4])):
        return None
    left = [int(closest_utr3[3]), cs]
    logging.debug('Left region: {}-{}'.format(left[0], left[1]))
    right = [cs, int(closest_utr3[4])]
    logging.debug('Right region: {}-{}'.format(right[0], right[1]))
    ccle_aln = pysam.AlignmentFile('/projects/btl/polya/sorted_star_ccle/cram_sorted_{}.cram'.format(best_site[0]), 'rc')
    ccle_exp_left = ccle_aln.count(chrom, left[0], left[1]) / (left[1] - left[0])
    logging.debug('Exp left: {}'.format(ccle_exp_left))
    ccle_exp_right = ccle_aln.count(chrom, right[0], right[1]) / (right[1] - right[0])
    logging.debug('Exp right: {}'.format(ccle_exp_right))
    ccle_diff = ccle_exp_left - ccle_exp_right
    logging.debug('CCLE diff (left - right): {}'.format(ccle_diff))
    for g in all_gtex:
        logging.debug('Comparing to GTEX {}'.format(g))
        #sprint('\r{}'.format(('\t').join([str(x) for x in g])))
        try:
            gtex_aln = pysam.AlignmentFile('/projects/btl/polya/gtex/star/crams/{}.star.cram'.format(g[2]), 'rc')
        except IOError:
            continue
        if not gtex_aln.has_index():
            logging.debug('indexing {} ...'.format(gtex_aln.filename))
            continue
            #pysam.index(gtex_aln.filename)
            #gtex_aln = pysam.AlignmentFile('/projects/btl/polya/gtex/star/crams/{}.star.cram'.format(g[2]), 'rc')
        gtex_exp_left = gtex_aln.count(chrom, left[0], left[1])
        logging.debug('Exp left: {}'.format(gtex_exp_left))
        gtex_exp_right = gtex_aln.count(chrom, right[0], right[1])
        logging.debug('Exp right: {}'.format(gtex_exp_right))
        gtex_diff = gtex_exp_left - gtex_exp_right
        logging.debug('GTEX diff (left - right): {}'.format(gtex_diff))
        #print '-'*20
        #print ccle_diff
        #print gtex_diff
        #print '-'*20
        if ((ccle_diff <= 0) and (gtex_diff >= 0)) or ((ccle_diff >= 0) and (gtex_diff <= 0)):
            #raw_input('Pausing...')
            logging.debug('Found difference in {}'.format(g))
            return g
    return None
        
# site format (string) [tab_delim]: 'CHROM GENE AVG_CLEAVAGE SCORE FREQ AVG_DIST'
def compare(site, gnovel):
    chrom = site[0]
    gene = site[1]
    if chrom not in gnovel:
        return None
    elif gene not in gnovel[chrom]:
        return None
    relevant = sorted(gnovel[chrom][gene], key=lambda x: int(x[17]) + int(x[18]))
    pos = best_site = pas = None
    for i, site in enumerate(relevant):
        try:
            pot_pas = int(site[15][1:-1].split(', ')[1])
            if (pas == None) or (pot_pas < pas):
                pas = pot_pas
                pos = i
                best_site = site[:]
        except IndexError:
            continue
    if not pas:
        best_site = relevant[0]
    return best_site

def pickGtex(bs, gtex, exclude=[]):
    tissue = bs[2].lower()
    picked = random.choice([x for x in gtex if x[0].lower() == tissue and x not in exclude])
    picked_colour = '{},{},{}'.format(random.randint(0,255),random.randint(0,255),random.randint(0,255))
    picked_url = 'http://bcgsc.ca/downloads/dmacmillan/{}.bw'.format(picked[2])
    picked_path = '/gsc/www/bcgsc.ca/downloads/dmacmillan/{}.bw'.format(picked[2])
    if not os.path.isfile(picked_path):
        shutil.copy('/projects/btl/polya/gtex_bigwigs/{}.bw'.format(picked[2]), '/gsc/www/bcgsc.ca/downloads/dmacmillan/')
    picked_track = 'track type=bigWig name={} description="bigwig track for {}" color="{}" visibility=full bigDataUrl={}\n'.format(picked[2], picked[2], picked_colour, picked_url)
    return picked, picked_track

def generateKleatTrack(ccle_id, strand):
    webdir = '/gsc/www/bcgsc.ca/downloads/dmacmillan/'
    if strand == '+':
        path = '/projects/dmacmillanprj2/polya/ccle/filteredKleats/kleats_plus_added/{}.plus.bg'.format(ccle_id)
    else:
        path = '/projects/dmacmillanprj2/polya/ccle/filteredKleats/kleats_plus_added/{}.minus.bg'.format(ccle_id)
    #if not os.path.isfile(os.path.join(webdir, os.path.basename(path))):    
        #sprint('\rCopying KLEAT track for {} ...'.format(ccle_id))
        #shutil.copy(path, webdir)
    f = open(path, 'r')
    data = ('').join(f.readlines()[1:])
    f.close()
    return 'track type=bedGraph name={}_novel_site description="Cleavage Sites for {} transcripts in {}" color="0,0,255" visibility=full\n{}\n'.format(ccle_id, strand, ccle_id, data)

def generateCcleTracks(*ccles):
    result = ''
    webdir = '/gsc/www/bcgsc.ca/downloads/dmacmillan/'
    for ccle in ccles:
        colour = '{},{},{}'.format(random.randint(0,255),random.randint(0,255),random.randint(0,255))
        ccle_bw_path = '/projects/btl/polya/ccle_bigwigs/{}.bw'.format(ccle)
        if not os.path.isfile('{}{}.bw'.format(webdir, ccle)):    
            #sprint('\rCopying bigwig track for {} ...'.format(ccle))
            shutil.copy(ccle_bw_path, webdir)
        big_data_url = 'http://bcgsc.ca/downloads/dmacmillan/{}.bw'.format(ccle)
        track_line = 'track type=bigWig name="{}" description="genomecov bigwig track for {}" color="{}" visibility=full bigDataUrl={}\n'.format(ccle, ccle, colour, big_data_url)
        result += track_line
    return result

def generateGtexTracks(*gtexs):
    result = ''
    webdir = '/gsc/www/bcgsc.ca/downloads/dmacmillan/'
    for gtex in gtexs:
        colour = '{},{},{}'.format(random.randint(0,255),random.randint(0,255),random.randint(0,255))
        gtex_bw_path = '/projects/btl/polya/gtex_bigwigs/{}.bw'.format(gtex)
        if not os.path.isfile('{}{}.bw'.format(webdir, ccle)):    
            #sprint('\rCopying bigwig track for {} ...'.format(gtex))
            shutil.copy(gtex_bw_path, webdir)
        big_data_url = 'http://bcgsc.ca/downloads/dmacmillan/{}.bw'.format(gtex)
        track_line = 'track type=bigWig name="{}" description="genomecov bigwig track for {}" color="{}" visibility=full bigDataUrl={}\n'.format(gtex, gtex, colour, big_data_url)
        result += track_line
    return result

def generateBrowserTrack(chrom, cs, window = 2000):
    start = cs - window
    end = cs + window
    browser = 'browser position {}:{}-{}\n'.format(chrom, start, end)
    browser += 'browser hide all\n'
    browser += 'browser pack knownGene\n'
    browser += 'browser pack refGene\n'
    browser += 'browser pack acembly\n'
    browser += 'browser pack ensGene\n'
    return browser

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

for clust in clusters:
    clust['best_site'] = None
    clust['best_score'] = 0
    logging.debug('Total sites in cluster: {}'.format(len(clust['sites'])))
    for i,site in enumerate(clust['sites']):
        logging.debug('site index: {}'.format(i))
        logging.debug('site: {}'.format(site))
#        try:
#            score = (int(site[16]) + int(site[17]) + int(site[18]) + int(site[19])) / int(site[15][1])
#        except (ValueError, IndexError) as e:
#            score = (int(site[16]) + int(site[17]) + int(site[18]) + int(site[19]))
#        clust['sites'][i].append(score)
#        if score > clust['best_score']:
#            clust['best_site'] = site
#            clust['best_score'] = score
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
        chrom = clust['best_site'][4]
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
        if diff <= args.threshold:
            logging.debug('Skipped because diff is too low')
            continue
        print '{}\t{}'.format(name, cs)
#    novel_track = 'track name="Novel Site" description="Novel site {}" color="255,50,50" visibility=full\n'.format(cs)
#    novel_track += '{}\t{}\t{}\t{}\n'.format(chrom, cs-1, cs, clust['best_score'])
#    regions_track = 'track name="UTR3 Regions" description="Regions of UTR3" visibility=full\n'
#    regions_track += '{}\t{}\t{}\t{}\n'.format(chrom, utr3_start, cs, left)
#    regions_track += '{}\t{}\t{}\t{}\n'.format(chrom, cs, utr3_end, right)
#    to_write = generateBrowserTrack(chrom, cs) + generateCcleTracks(name) + novel_track + regions_track
#    fname = '{}_{}'.format(name, cs)
#    with open(os.path.join(webdir, fname), 'w') as f:
#        f.write(to_write)
#    url = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hgt.customText=http://bcgsc.ca/downloads/dmacmillan/{}&hgt.psOutput=on'.format(fname)
#    page = requests.get(url)
#    soup = BeautifulSoup(page.content, 'html.parser')
#    pdf = soup.find_all('a')[-7]['href']
#    pdf = 'http://genome.ucsc.edu/cgi-bin/' + pdf
#    urllib.urlretrieve(pdf, os.path.join(args.outdir, fname + '.pdf'))
