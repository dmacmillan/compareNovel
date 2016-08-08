import argparse, os, sys, logging
import random
import shutil
import pysam

parser = argparse.ArgumentParser(description='')

parser.add_argument('all_novel', help='A file containing information for all the novel sites')
parser.add_argument('strong_novel', help='A file containing novel sites and their relative frequencies')
parser.add_argument('gtex_list', help='File containing information on selected GTEX samples')
parser.add_argument('kleats_filtered', nargs='+', help='Filtered KLEAT files')
parser.add_argument('-a', '--annotation', default='/projects/dmacmillanprj2/polya/ccle/novel_cleavage_events/pysam_version/with_utr3/ucsc.utr3.gtf', help='UCSC-downloaded annotation file')
parser.add_argument("-l", "--log", dest="logLevel", default='WARNING', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help='Set the logging level. Default = "WARNING"')
parser.add_argument('-n', '--name', default='result', help='Name for the file output. Default is "result"')
parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is {}'.format(os.getcwd()))

args = parser.parse_args()

filtered_kleats_path = '/projects/dmacmillanprj2/polya/ccle/filteredKleats/kleats_plus_added'
webdir = '/gsc/www/bcgsc.ca/downloads/dmacmillan/'

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
    novel = [x.strip().split('\t') for x in f.readlines()]
    novel = novel[1:]

gnovel = {}

for i in novel:
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
    chrom = best_site[4]
    closest_utr3 = getClosestUtr3(best_site, utrs)
    #print closest_utr3
    cs = int(best_site[6])
    if not (int(closest_utr3[3]) < cs < int(closest_utr3[4])):
        return None
    left = [int(closest_utr3[3]), cs]
    right = [cs, int(closest_utr3[4])]
    ccle_aln = pysam.AlignmentFile('/projects/btl/polya/sorted_star_ccle/cram_sorted_{}.cram'.format(best_site[0]), 'rc')
    ccle_exp_left = ccle_aln.count(chrom, left[0], left[1])
    ccle_exp_right = ccle_aln.count(chrom, right[0], right[1])
    ccle_diff = ccle_exp_left - ccle_exp_right
    logging.debug('Looking through gtex ...')
    for g in all_gtex:
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
        gtex_exp_right = gtex_aln.count(chrom, right[0], right[1])
        gtex_diff = gtex_exp_left - gtex_exp_right
        #print '-'*20
        #print ccle_diff
        #print gtex_diff
        #print '-'*20
        if ((ccle_diff <= 0) and (gtex_diff >= 0)) or ((ccle_diff >= 0) and (gtex_diff <= 0)):
            raw_input('Pausing...')
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

def getTracks(bs, all_gtex, window=2000):
    try:
        pas = [int(x) for x in bs[15][1:-1].split(', ')]
    except IndexError as e:
        pas = None
    chrom = bs[4]
    start = int(bs[6]) - window
    end = int(bs[6]) + window
    strength = int(bs[16]) + int(bs[17]) + int(bs[18]) + int(bs[19])
    name = bs[0]
    strand = bs[10]
    if strand == '-':
        kleat_track = '/projects/btl/polya/ccle/{}/kleat/{}.minus.bg'.format(name, name)
    else:
        kleat_track = '/projects/btl/polya/ccle/{}/kleat/{}.plus.bg'.format(name, name)
    kleat_track = filterKleatTrack(kleat_track)
    colour = '{},{},{}'.format(random.randint(0,255),random.randint(0,255),random.randint(0,255))
    description='{}_bigwig'.format(name)
    readman_loc='/projects/btl/polya/ccle_bigwigs/{}.bw'.format(name)
    bigDataUrl='http://bcgsc.ca/downloads/dmacmillan/{}.bw'.format(name)
    browser = 'browser position {}:{}-{}\n'.format(chrom, start, end)
    browser += 'browser hide all\n'
    browser += 'browser pack knownGene\n'
    browser += 'browser pack refGene\n'
    browser += 'browser pack acembly\n'
    browser += 'browser pack ensGene\n'
    bigwig_track = 'track type=bigWig name={} description="{}" color="{}" visibility=full bigDataUrl={}\n'.format(name, description, colour, bigDataUrl)
    novel_track = 'track type=bedGraph name={}_novel_site description="Novel cleavage event at {}" color="0,0,255" visibility=full\n'.format(name, bs[6])
    novel_track += ('\t').join([chrom, bs[6], str(int(bs[6])+1), str(strength)]) + '\n'
    gtex, gtex_track = pickGtex(bs, all_gtex)
    if not os.path.isfile('/gsc/www/bcgsc.ca/downloads/dmacmillan/{}.bw'.format(name)):    
        shutil.copy(readman_loc, '/gsc/www/bcgsc.ca/downloads/dmacmillan/')
    output = '/gsc/www/bcgsc.ca/downloads/dmacmillan/{}_{}'.format(name, bs[6])
    with open(output, 'w') as f:
        f.write(browser + bigwig_track + novel_track + gtex_track)
        if pas:
            pas_track = 'track name={}_pas description="PAS with Strength {} and position {}" color="255,0,0" visibility=full\n'.format(name, pas[1], pas[0])
            pas_track += '{}\t{}\t{}\n'.format(chrom, pas[0], pas[0]+6)
            f.write(pas_track)
    #return output
    return 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hgt.customText=http://bcgsc.ca/downloads/dmacmillan/{}_{}'.format(name, bs[6])

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

def filterKleatTrack(kleat_track):
    return None

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

#for x in strong:
#    getTracks(compare(x, gnovel), gtex)
#print getTracks(compare(strong[0], gnovel), all_gtex)

for st in strong[2:]:
    res = detectSwitch(st, all_gtex, utrs, gnovel)
    if res:
        best_match = compare(st, gnovel)
        try:
            pas = [int(x) for x in best_match[15][1:-1].split(', ')]
        except (ValueError, IndexError) as e:
            pas = None
        novel_track = 'track name="Novel Site" description="Novel site {}" color="255,50,50" visibility=full\n'.format(st[2])
        novel_track += '{}\t{}\t{}\t{}\n'.format(st[0], st[2], int(st[2])+1, st[3])
        ccle = generateCcleTracks(best_match[0])
        gtex = generateGtexTracks(res[2])
        kleat = generateKleatTrack(best_match[0], best_match[10])
        browser = generateBrowserTrack(best_match[4], int(best_match[6]))
        pas_track = ''
        if pas:
            pas_track = 'track name={}_pas description="{}" color="255,0,0" visibility=full\n'.format(best_match[0], pas)
            pas_track += '{}\t{}\t{}\n'.format(st[0], pas[0], pas[0]+6)
        to_write = browser + ccle + pas_track + kleat + novel_track + gtex
        fname = '{}_{}_{}'.format(st[0], st[1], st[2])
        with open(os.path.join(webdir, fname), 'w') as f:
            f.write(to_write)
        print 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hgt.customText=http://bcgsc.ca/downloads/dmacmillan/{}'.format(fname)

    #if res:
    #    getTracks(compare(x, gnovel), all_gtex)
#detectSwitch(strong[0], all_gtex, utrs, gnovel)
