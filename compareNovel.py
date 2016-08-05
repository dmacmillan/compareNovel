import argparse, os, sys, logging
import random
import subprocess
import pysam

parser = argparse.ArgumentParser(description='')

parser.add_argument('all_novel', help='A file containing information for all the novel sites')
parser.add_argument('strong_novel', help='A file containing novel sites and their relative frequencies')
parser.add_argument('gtex_list', help='File containing information on selected GTEX samples')
parser.add_argument('-a', '--annotation', default='/projects/dmacmillanprj2/polya/ccle/novel_cleavage_events/pysam_version/with_utr3/ucsc.utr3.gtf', help='UCSC-downloaded annotation file')
parser.add_argument("-l", "--log", dest="logLevel", default='WARNING', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help='Set the logging level. Default = "WARNING"')
parser.add_argument('-n', '--name', default='result', help='Name for the file output. Default is "result"')
parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is {}'.format(os.getcwd()))

args = parser.parse_args()

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
            clostest_utr3 = utrs[i][:]
    return closest_utr3

def detectSwitch(site, gtex, utrs):
    closest_utr3 = getClosestUtr3(site, utrs)
    cs = int(site[6])
    left = [int(closest_utr3[3]), cs]
    right = [cs, int(closest_utr3[4])]
    bw = pysam.AlignmentFile('/projects/btl/polya/sorted_star_ccle/cram_sorted_{}.cram'.format(site[]), 'rc')

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

def getUrl(best_site, db='hg19', track_type='bigWig', visibility='full', window=2000):
    chrom = best_site[4]
    start = int(best_site[6]) - window
    end = int(best_site[6]) + window
    colour = '%22{},{},{}%22'.format(random.randint(0,255),random.randint(0,255),random.randint(0,255))
    description='{}_bigwig'.format(best_site[0])
    readman_loc='/projects/btl/polya/ccle_bigwigs/{}.bw'.format(best_site[0])
    bigDataUrl='/gsc/www/bcgsc.ca/downloads/dmacmillan/{}.bw'.format(best_site[0])
    # whitespace == %20
    # quotes == %22
    if not os.path.isfile('/gsc/www/bcgsc.ca/downloads/dmacmillan/{}.bw'.format(best_site[0])):    
        subprocess.call(['cp', readman_loc, '/gsc/www/bcgsc.ca/downloads/dmacmillan/'])
    return 'http://genome.ucsc.edu/cgi-bin/hgTracks?db={}&position={}:{}-{}&knownGene=pack&refGene=pack&acembly=pack&ensGene=pack&{}=full&hgct_customText=track%20type={}%20name={}%20description=%22{}%22%20visibility={}%20color={}%20bigDataUrl={}'.format(db, chrom, start, end, best_site[0], track_type, best_site[0], description, visibility, colour, bigDataUrl)

def pickGtex(bs, gtex, exclude=[]):
    tissue = bs[2].lower()
    picked = random.choice([x for x in gtex if x[0].lower() == tissue and x not in exclude])
    picked_colour = '{},{},{}'.format(random.randint(0,255),random.randint(0,255),random.randint(0,255))
    picked_url = 'http://bcgsc.ca/downloads/dmacmillan/{}.bw'.format(picked[2])
    picked_path = '/gsc/www/bcgsc.ca/downloads/dmacmillan/{}.bw'.format(picked[2])
    if not os.path.isfile(picked_path):
        subprocess.call(['cp', '/projects/btl/polya/gtex_bigwigs/{}.bw'.format(picked[2]), '/gsc/www/bcgsc.ca/downloads/dmacmillan/'])
    picked_track = 'track type=bigWig name={} description="bigwig track for {}" color="{}" visibility=full bigDataUrl={}\n'.format(picked[2], picked[2], picked_colour, picked_url)
    return picked, picked_track

def getTracks(bs, all_gtex, window=2000):
    try:
        pas = [int(x) for x in bs[15][1:-1].split(', ')]
    except IndexError:
        pas = None
    chrom = bs[4]
    start = int(bs[6]) - window
    end = int(bs[6]) + window
    strength = int(bs[16]) + int(bs[17]) + int(bs[18]) + int(bs[19])
    name = bs[0]
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
        subprocess.call(['cp', readman_loc, '/gsc/www/bcgsc.ca/downloads/dmacmillan/'])
    output = '/gsc/www/bcgsc.ca/downloads/dmacmillan/{}_{}'.format(name, bs[6])
    with open(output, 'w') as f:
        f.write(browser + bigwig_track + novel_track + gtex_track)
        if pas:
            pas_track = 'track name={}_pas description="PAS with Strength {} and position {}" color="255,0,0" visibility=full\n'.format(name, pas[1], pas[0])
            pas_track += '{}\t{}\t{}\n'.format(chrom, pas[0], pas[0]+6)
            f.write(pas_track)
    #return output
    return 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hgt.customText=http://bcgsc.ca/downloads/dmacmillan/{}_{}'.format(name, bs[6])

#for x in strong:
#    getTracks(compare(x, gnovel), gtex)
print getTracks(compare(strong[0], gnovel), all_gtex)
