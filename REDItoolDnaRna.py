#!/home/epicardi/bin/python27/bin/python
# Copyright (c) 2013-2014 Ernesto Picardi <ernesto.picardi@uniba.it>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import sys, os, time, math, random, getopt, operator, string, errno
try: import pysam
except: sys.exit('Pysam module not found.')
from multiprocessing import Process, Queue
from Queue import Empty
import gzip

version='1.0'
pid=str(os.getpid()+random.randint(0,999999999))

def usage():
	print """
USAGE: python REDItoolDnaRNA.py [options]
Options:
-i		RNA-Seq BAM file
-j		DNA-Seq BAM file(s separated by comma) or folder
-I		Sort input RNA-Seq BAM file
-J		Sort input DNA-Seq BAM file
-f		Reference in fasta file
-C		Base interval to explore [100000]
-k		List of chromosomes to skip separated by comma or file
-t		Number of threads [1]
-o		Output folder [rediFolder_%s]
-F		Internal folder name [null]
-M		Save a list of columns with quality scores
-c		Min. read coverage (dna,rna) [10,10]
-Q		Fastq offset value (dna,rna) [33,33]
-q		Min. quality score (dna,rna) [25,25]
-m		Min. mapping quality score (dna,rna) [25,25]
-O		Min. homoplymeric length (dna,rna) [5,5]
-s		Infer strand (for strand oriented reads) [1]
-g		Strand inference type 1:maxValue 2:useConfidence [1]
-x		Strand confidence [0.70]
-S		Strand correction
-G		Infer strand by GFF annotation (must be GFF and sorted, otherwise use -X)
-K		GFF File with positions to exclude (must be GFF and sorted, otherwise use -X)
-T		Work only on given GFF positions (must be GFF and sorted, otherwise use -X)
-X		Sort annotation files
-e		Exclude multi hits in RNA-Seq
-E		Exclude multi hits in DNA-Seq
-d		Exclude duplicates in RNA-Seq
-D		Exclude duplicates in DNA-Seq
-p		Use paired concardant reads only in RNA-Seq
-P		Use paired concardant reads only in DNA-Seq
-u		Consider mapping quality in RNA-Seq
-U		Consider mapping quality in DNA-Seq
-a		Trim x bases up and y bases down per read [0-0] in RNA-Seq
-A		Trim x bases up and y bases down per read [0-0] in DNA-Seq
-b		Blat folder for correction in RNA-Seq
-B		Blat folder for correction in DNA-Seq
-l		Remove substitutions in homopolymeric regions in RNA-Seq
-L		Remove substitutions in homopolymeric regions in DNA-Seq
-v		Min. num. of reads supporting the variation [3] for RNA-Seq
-n		Min. editing frequency [0.1] for RNA-Seq
-N		Min. variation frequency [0.1] for DNA-Seq
-z		Exclude positions with multiple changes in RNA-Seq
-Z		Exclude positions with multiple changes in DNA-Seq
-W		Select RNA-Seq positions with defined changes (separated by comma ex: AG,TC) [default all]
-R		Exclude invariant RNA-Seq positions
-V		Exclude sites not supported by DNA-Seq
-w		File containing splice sites annotations
-r		Num. of bases near splice sites to explore [4]
--gzip	Gzip output files
-h		Print this help
--help

"""%(pid)

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:f:k:t:o:c:Q:q:m:O:s:edpuA:a:B:b:lLv:n:EPr:hIXG:K:j:C:JDUzw:N:ZW:RVMT:F:x:g:S",["help","gzip"])
except getopt.GetoptError as err:
	print str(err) # will print something like "option -a not recognized"
	usage()
	sys.exit(2)

strconf=0.70 #confidenza strand
useconf=0
corrstr=0
bamfile=''
gbamfile=[]
dgbamfile={} # dizionario di gbamfile
fastafile=''
sortbam=0
sortgbam=0
nochrs=[]
NCPU=1
infolder=''
outfolder_='rediFolder_%s' %(pid)
MINCOV=10
QVAL=33
MQUAL=25
MAPQ=20
homo=5
rmpv = '0-0'
rmp = [int(x) for x in rmpv.split('-')]
getstrand=0 # considera la strand
exh=0 # escludi multi hits
exd=0 # escludi duplicati
conc=0 # se presenti paired-end, usa solo quelle concordanti
mq=0 # considera il map quality
rmnuc=0 # rimuovi nucleotide a monte ed a valle delle read; connesso a rmp e rmpv
blatr=0 # applica la correzione blat
blatfolder=''
rmsh=0 # rimuovi sostituzioni in omopolimeri di lunghezza maggiore o uguale a homo
vnuc=3 # numero minimo di basi che supportano la variazione
mmf=0.1 # frequenza minima della variazione
exms=0 # escludi sostituzioni multiple
exss=0 # escludi posizioni introniche nei pressi dei siti di splicing a nss nucleotidi 
nss=4 # basi introniche da esplorare per ogni sito si splicing
splicefile='' #'splicesites.hg18.sorted.txt'
usubs=[x+y for x in 'ACGT' for y in 'ACGT' if x!=y] # use these substitutions [default all]
annfile='' # use annotation file for strand correction and features
sortann=0 # sort annotation file
uann=0 # use annotation
exfile='' # use annotations to exclude positions
expos=0 #
chunckval=100000
exinv=0
slist=0
slistfile=''
gslistfile=''
wfile='' # working file. GFF annotations to use
uwf=0 # use working file
###########################
# for DNA-Seq
###########################
gMINCOV=10
gQVAL=33
gMQUAL=25
gMAPQ=20
ghomo=5
grmpv = '0-0'
grmp = [int(x) for x in rmpv.split('-')]
gexh=0 # escludi multi hits
gexd=0 # escludi duplicati
gconc=0 # se presenti paired-end, usa solo quelle concordanti
gmq=0 # considera il map quality
grmnuc=0 # rimuovi nucleotide a monte ed a valle delle read; connesso a rmp e rmpv
gblatr=0 # applica la correzione blat
gblatfolder=''
grmsh=0 # rimuovi sostituzioni in omopolimeri di lunghezza maggiore o uguale a ghomo
gmmf=0.1 # frequenza minima della variazione
exnonh=0 # escludi posizioni non omozigoti
exnosupp=0 # escludi posizioni non supportate da DNA-Seq
nogbam=0 # esiste il bam per DNA-Seq
unchange1=1
unchange2=0
gziptab=0

for o, a in opts:
	if o in ("-h","--help"):
		usage()
		sys.exit()
	elif o == "-i": bamfile=a
	elif o == "-j":
		if os.path.isdir(a): gbamfile=[(os.path.join(a,x),0) for x in os.listdir(a) if x[-4:]=='.bam']
		else: gbamfile=[(x,0) for x in a.split(',') if x.strip()!='']
		dgbamfile=dict(gbamfile)
	elif o == "-f": fastafile=a
	elif o == "-k":
		if os.path.exists(a):
			f=open(a)
			nochrs=[x.strip() for x in f if x.strip()!='']
			f.close()
		else: nochrs = a.split(',')
	elif o == "-t": NCPU=int(a)
	elif o == "-F": infolder=a	
	elif o == "-o": outfolder_=a
	elif o == "-c":
		MINCOV=int(a.split(',')[1])
		gMINCOV=int(a.split(',')[0])
	elif o == "-Q":
		QVAL=int(a.split(',')[1])
		gQVAL=int(a.split(',')[0])
	elif o == "-q":
		MQUAL=int(a.split(',')[1])
		gMQUAL=int(a.split(',')[0])
	elif o == "-m":
		MAPQ=int(a.split(',')[1])
		gMAPQ=int(a.split(',')[0])
	elif o == "-O":
		homo=int(a.split(',')[1])
		ghomo=int(a.split(',')[0])
	elif o == "-S": corrstr=1
	elif o == "-x": strconf=float(a)
	elif o == "-g":
		if a=='2': useconf=1
	elif o == "-s":
		getstrand=1
		if int(a)==1: unchange1,unchange2=1,0
		elif int(a)==0: unchange1,unchange2=0,0
		elif int(a)==2: unchange1,unchange2=0,1
		elif int(a)==12: unchange1,unchange2=1,1
	elif o == "-e": exh=1
	elif o == "-V": exnosupp=1
	elif o == "-E": gexh=1
	elif o == "-d": exd=1
	elif o == "-D": gexd=1
	elif o == "-p": conc=1
	elif o == "-P": gconc=1
	elif o == "-I": sortbam=1
	elif o == "-J": sortgbam=1
	elif o == "-X": sortann=1
	elif o == "-R": exinv=1
	elif o == "-C": chunckval=int(a)
	elif o == "-u": mq=1
	elif o == "-U": gmq=1
	elif o == "-M": slist=1
	elif o == "-a":
		rmpv = a
		try:
			rmp = [int(x) for x in rmpv.split('-')]
			rmnuc=1
		except: rmnuc=0
	elif o == "-A":
		grmpv = a
		try:
			grmp = [int(x) for x in grmpv.split('-')]
			grmnuc=1
		except: grmnuc=0
	elif o == "-b":
		blatfolder=a
		if os.path.exists(blatfolder): blatr=1
	elif o == "-B":
		gblatfolder=a
		if os.path.exists(gblatfolder): gblatr=1
	elif o == "-l": rmsh=1
	elif o == "-L": grmsh=1
	elif o == "-v": vnuc=int(a)
	elif o == "-n": mmf=float(a)
	elif o == "-N": gmmf=float(a)
	elif o == "-z": exms=1
	elif o == "-Z": exnonh=1
	elif o == "-W": usubs=[x.upper() for x in a.split(',') if x.strip()!='']
	elif o == "-w":
		splicefile=a
		if os.path.exists(splicefile): exss=1
	elif o == "-K":
		exfile=a
		if os.path.exists(exfile): expos=1
	elif o == "-T":
		wfile=a
		if os.path.exists(wfile): uwf=1	
	elif o == "-r": nss=int(a)
	elif o == "-G":
		annfile=a
		uann=1
	elif o == "--gzip": gziptab=1
	else:
		assert False, "Unhandled Option"

#######
commandLine=' '.join(sys.argv[1:])
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
params=[]
#Input parameters
params.append('REDItoolDnaRna version %s\n' %(version))
params.append('User command line: %s\n' %(commandLine))
params.append('Analysis ID: %s\n' %(pid))
params.append('Analysis time: %s\n' %(script_time))
params.append('-i --> RNA-Seq BAM file: %s\n' %(bamfile))
params.append('-j --> DNA-Seq BAM file(s): %s\n' %(','.join(dgbamfile.keys())))
params.append('-I --> Sort RNA-Seq BAM file: %i\n' %(sortbam))
params.append('-J --> Sort DNA-Seq BAM file: %i\n' %(sortgbam))
params.append('-f --> Reference file: %s\n' %(fastafile))
params.append('-C --> Base interval to explore: %i\n' %(chunckval))
params.append('-k --> Regions to exclude: %s\n' %(','.join(nochrs)))
params.append('-t --> Number of working threads: %i\n' %(NCPU))
params.append('-o --> Output folder: %s\n' %(outfolder_))
params.append('-F --> Infolder folder: %s\n' %(infolder))
params.append('-M --> Save a list of columns with quality scores: %i\n' %(slist))
params.append('-c --> Min. per base coverage DNA-RNA: %i-%i\n' %(MINCOV,gMINCOV))
params.append('-Q --> FastQ offset value DNA-RNA: %i-%i\n' %(gQVAL,QVAL))
params.append('-q --> Min. per base quality DNA-RNA: %i-%i\n' %(gMQUAL,MQUAL))
params.append('-m --> Min. mapping quality DNA-RNA: %i-%i\n' %(gMAPQ,MAPQ))
params.append('-O --> Min. homoplymeric length DNA-RNA: %i-%i\n' %(ghomo,homo))
params.append('-s --> Infer strand: %i - %i-%i\n' %(getstrand,unchange1,unchange2))
params.append('-g --> Use confidence: %i\n' %(useconf))
params.append('-x --> Strand confidence: %.2f\n' %(strconf))
params.append('-S --> Strand correction : %i\n' %(corrstr))
params.append('-G --> GFF annotation to infer strand: %s\n' %(annfile))
params.append('-K --> File with positions to exclude: %s\n' %(exfile))
params.append('-T --> Work only on given GFF positions: %s\n' %(wfile))
params.append('-X --> Sort annotation files: %i\n' %(sortann))
params.append('-e --> Exclude multi hits in RNA-Seq: %i\n' %(exh))
params.append('-E --> Exclude multi hits in DNA-Seq: %i\n' %(gexh))
params.append('-d --> Exclude duplicates in RNA-Seq: %i\n' %(exd))
params.append('-D --> Exclude duplicates in DNA-Seq: %i\n' %(gexd))
params.append('-p --> Use paired concardant reads in RNA-Seq: %i\n' %(conc))
params.append('-P --> Use paired concardant reads in DNA-Seq: %i\n' %(gconc))
params.append('-u --> Consider mapping quality in RNA-Seq: %i\n' %(mq))
params.append('-U --> Consider mapping quality in DNA-Seq: %i\n' %(gmq))
params.append('-a --> Trim x bases up and y bases down per RNA read: %i - %i-%i\n' %(rmnuc,rmp[0],rmp[1]))
params.append('-A --> Trim x bases up and y bases down per DNA read: %i - %i-%i\n' %(grmnuc,grmp[0],grmp[1]))
params.append('-b --> Blat folder for correction in RNA-Seq: %s\n' %(blatfolder))
params.append('-B --> Blat folder for correction in DNA-Seq: %s\n' %(gblatfolder))
params.append('-l --> Remove substitutions in homopolymeric regions for RNA-Seq: %i\n' %(rmsh))
params.append('-L --> Remove substitutions in homopolymeric regions for DNA-Seq: %i\n' %(grmsh))
params.append('-v --> Min. num. of reads supporting the variation: %i\n' %(vnuc))
params.append('-n --> Min. editing frequency for RNA-Seq: %.2f\n' %(mmf))
params.append('-N --> Min. editing frequency for DNA-Seq: %.2f\n' %(gmmf))
params.append('-z --> Exclude positions with multiple changes in RNA-Seq: %i\n' %(exms))
params.append('-Z --> Exclude positions with multiple changes in DNA-Seq: %i\n' %(exnonh))
params.append('-W --> Select RNA-Seq positions with defined changes: %s\n' %(','.join(usubs)))
params.append('-R --> Exclude invariant RNA-Seq positions: %i\n' %(exinv))
params.append('-V --> Exclude sites not supported by DNA-Seq: %i\n' %(exnosupp))
params.append('-w --> File containing splice sites annotations: %s\n' %(splicefile))
params.append('-r --> Num. of bases near splice sites to explore: %i\n' %(nss))
params.append('--gzip --> Gzip output files: %i\n' %(gziptab))
#######
        
def get_no(pvalue,siglevel,ngenes): # No Correction
	lista=[]
	pp=siglevel
	y=0
	for i in pvalue:
		p=i[0]
		if p<=siglevel:
			lista.append(i)
			y+=1
	return lista,y,pp

def get_b(pvalue,siglevel,ngenes): # Bonferroni
	pvalue.sort()
	lista=[]
	y=0
	#bcorr=siglevel/ngenes
	pp=1.0
	for i in pvalue:
		p=i[0]*ngenes
		if p<=siglevel:
			lista.append(i)
			#lista[i[1]]=i[0]
			y+=1
			if p<pp: pp=p
	#print "Passed:",y,pp
	return lista,y,pp

def get_bh(pvalue,siglevel,ngenes): # B-H
	pvalue.sort()
	#print ngenes
	lista=[]
	x=1
	y=0
	p=0
	for i in pvalue:
		nf=i[0]*ngenes
		fdr=nf/x
		if fdr<=siglevel:
			#dic[i[1]]=i[0]
			lista.append(i)
			p=i[0]
			y+=1
		x+=1
	#print "Passed:",y,p
	return lista,y,p

def getTail(pp):
	if ftail=='l': return pp.left_tail
	elif ftail=='r': return pp.right_tail
	elif ftail=='t': return pp.two_tail
	
def getDicSS(dicp): # dicp = dizionario con le frequenze di sostituzione
	dicpp={}
	for i in dicp:
		if i[0]!=i[1]:
			dicpp[i]=1-dicp[i]
	return dicpp

def getFreads(bases):
	fread={'A':0,'C':0,'G':0,'T':0}
	for i in range(4):
		if i==0: fread['A']=bases[i]
		elif i==1: fread['C']=bases[i]
		elif i==2: fread['G']=bases[i]
		elif i==3: fread['T']=bases[i]
	return fread

def getSub(ref,fread,dics):
	#fread={A,C,G,T}
	nref=fread[ref.upper()]
	sub=[(ref.upper()+i,nref,fread[i]) for i in fread if i!=ref.upper() and fread[i]!=0]
	allsub=' '.join([x[0] for x in sub])
	# lista del tipo [('AT', 50, 10), ('AG', 50, 2)]
	res=[] #[(int(dics[i[0]]*(i[1]+i[2])),((i[1]+i[2])-exp1),pvalue(i[1],i[2],int(dics[i[0]]*(i[1]+i[2])),((i[1]+i[2])-exp1))) for i in sub]
	for i in sub:
		#if binomial:
		#	pval=bdtrc(i[2],i[1]+i[2],(1.0-dics[i[0]]))
		#	#pval=Bprob(i[2],i[1]+i[2],(1.0-dics[i[0]]))
		#	#print i[2],i[1]+i[2],(1.0-dics[i[0]]),pval
		#	obs1,obs2,exp1,exp2=0,0,0,0
		obs1=i[1]
		obs2=i[2]
		exp1=int(dics[i[0]]*(i[1]+i[2]))
		exp2=((i[1]+i[2]) - exp1)
		pval=pvalue(obs1,obs2,exp1,exp2)
		pval=getTail(pval)
		res.append((i[0],obs1,obs2,exp1,exp2,str(pval)))
	if len(res)==1: return res[0][5] #,allsub,fread
	elif len(res) > 1:
		rr=[float(x[-1]) for x in res]
		idx=rr.index(min(rr))
		return res[idx][5] #,allsub,fread
	else: return '1.0' #,0,0
	
def BaseCount(seq,ref,mfr,VNUC):
	b={'A':0,'C':0,'G':0,'T':0}
	subs=[]
	subv=[]
	for i in seq.upper():
		if b.has_key(i): b[i]+=1
	for i in b:
		if not b.has_key(ref): continue
		if b[i]!=0 and i!=ref:
			vv=float(b[i])/(b[i]+b[ref])
			subv.append((b[i],vv,ref+i))
	subv.sort()
	subv.reverse()
	for i in subv:
		if i[0]>=VNUC and i[1]>=mfr: subs.append(i[2])
	freq=0.0
	if len(subs)==0: subs.append('-')
	else: freq=subv[0][1]	
	return sum(b.values()),[b['A'],b['C'],b['G'],b['T']],' '.join(subs),'%.2f'%(freq)

def meanq(v,n):
	try:m=float(v)/n
	except: m=0.0
	return '%.2f'%(m)
	
def rmHomo(sequp,seqdw,gh,ref):
	if len(sequp)==0 and len(seqdw)==0: return 0
	up,dw=0,0
	for i in seqdw:
		if i==ref:dw+=1
		else:break
	for i in sequp[::-1]:
		if i==ref:up+=1
		else:break
	hlen=up+dw+1
	if hlen >= gh : return 1
	return 0

def prop(tot,va):
	try: av=float(va)/tot
	except: av=0.0
	return av

def vstand(strand):
	vv=[(strand.count('+'),'+'),(strand.count('-'),'-'),(strand.count('*'),'*')]
	if vv[0][0]==0 and vv[1][0]==0: return '*'
	if useconf:
		totvv=sum([x[0] for x in vv[:2]])
		if prop(totvv,vv[0][0])>=strconf: return '+'
		if prop(totvv,vv[1][0])>=strconf: return '-'
		return '*'
	else:
		if vv[0][0]==vv[1][0] and vv[2][0]==0: return '+'
		return max(vv)[1]

def comp(s):
	a={'A':'T','T':'A','C':'G','G':'C'}
	ss=''
	for i in s.upper():
		if a.has_key(i): ss+=a[i]
		elif i==' ': ss+=' '
		elif i=='-': ss+='-'
		else: ss+='N'
	return ss	

def comp2(s):
	ss={}
	a={'A':'T','T':'A','C':'G','G':'C'}
	for i,j in enumerate('ACGT'): ss[a[j]]=s[i]
	return str([ss['A'],ss['C'],ss['G'],ss['T']])

def whereis(program):
	for path in os.environ.get('PATH', '').split(':'):
		if os.path.exists(os.path.join(path, program)) and not os.path.isdir(os.path.join(path, program)): return 1
	return 0

def vstrand(lista):
	if len(lista)==0: return '2'
	p=lista.count('+')
	m=lista.count('-')
	if p==len(lista): return '1'
	elif m==len(lista): return '0'
	else: return '2'

def getd(lines):
	d={}
	for i in lines:
		l=i.split('\t')
		if len(l)>=3:
			if l[2]=='+': strand='1'
			elif l[2]=='-': strand='0'
			else: strand='2'
		else: strand='2'
		d[int(l[1])]=strand
	return d

def checkSubs(s):
	if s=='-': return 1
	for i in s.split():
		if i in usubs: return 1
	return 0

def makeCluster(allcoord):
	cluster=[]
	remaining=[]
	c1=allcoord[0][0]
	c2=allcoord[0][1]
	for i in range(len(allcoord)):
		if allcoord[i]!=(c1,c2):
			if c1<=allcoord[i][0]<=c2:
				cluster.append(allcoord[i])
				if allcoord[i][1]>c2:
					c2=allcoord[i][1]
			else:
				remaining.append(allcoord[i])
		else:
			cluster.append((c1,c2))
	return (c1,c2),remaining

def newCoords(interval,start,end):
	coords=[]
	interval.sort()
	while len(interval)!=0:
		coord,interval=makeCluster(interval)
		coords.append(coord)
	c1,c2=coords[0][0],coords[-1][1]
	if c1 < start: c1=start
	if c2>end: c2=end
	return coords,c1,c2

def checkPos(coords,pos):
	for i in coords:
		if i[0]<=pos<=i[1]: return 1
	return 0

def parseFeat(line):
	l=line.split('\t')
	cc=(int(l[3])-1,int(l[4])-1)
	return cc

def normByStrand(seq_,strand_,squal_,mystrand_):
	st='+'
	if mystrand_=='0': st='-'
	seq,qual,squal='',0,''
	for i in range(len(seq_)):
		if strand_[i]==st:
			seq+=seq_[i]
			qual+=ord(squal_[i])-QVAL
			squal+=squal_[i]
	return seq,qual,squal

def normByBlat(seq_,strand_,squal_,blatc_,qqval):
	seq,qual,squal,strand='',0,'',''
	for i in range(len(seq_)):
		if blatc_[i]=='1':
			seq+=seq_[i]
			qual+=ord(squal_[i])-qqval
			squal+=squal_[i]
			strand+=strand_[i]
	return seq,qual,squal,strand
	
def testBlat(blc):
	if blc.count('1') > blc.count('0'): return 1
	return 0

###########################################################
###########################################################
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stderr.write("Script time --> START: %s\n"%(script_time))
sys.stderr.write("Analysis ID: %s\n"%(pid))
###########################################################

if not os.path.exists(bamfile):
	usage()
	sys.exit('RNA-Seq BAM file %s not found.' %(bamfile))
if sortbam:
	sys.stderr.write('Sorting RNA-Seq BAM file.\n')
	pysam.sort(bamfile,'sorted_%s'%(pid))
	os.rename(bamfile,bamfile+'_old')
	os.rename('sorted_%s.bam'%(pid),bamfile)
	sys.stderr.write('Indexing RNA-Seq BAM file.\n')
	pysam.index(bamfile)
if not os.path.exists(bamfile+'.bai') and not sortbam:
	sys.stderr.write('Indexing RNA-Seq BAM file.\n')
	pysam.index(bamfile)
###########################################################
dgdic={} # dizionario chr:bam file
for i in gbamfile:
	if not os.path.exists(i[0]):
		sys.stderr.write('DNA-Seq BAM file %s not found.\n' %(i[0]))
		sys.stderr.write('Working without DNA-Seq BAM file %s.\n' %(i[0]))
		del dgbamfile[i[0]]
	else:
		if sortgbam:
			sys.stderr.write('Sorting DNA-Seq BAM file %s.\n' %(i[0]))
			pysam.sort(i[0],'sorted_%s'%(pid))
			os.rename(i[0],i[0]+'_old')
			os.rename('sorted_%s.bam'%(pid),i[0])
			sys.stderr.write('Indexing DNA-Seq BAM file %s.\n' %(i[0]))
			pysam.index(i[0])
		if not os.path.exists(i[0]+'.bai') and not sortgbam:
			sys.stderr.write('Indexing DNA-Seq BAM file %s.\n' %(i[0]))
			pysam.index(i[0])
if len(gbamfile)==0:
	sys.stderr.write('Working without DNA-Seq BAM file(s).\n')
	nogbam=1
else:
	for i in dgbamfile:
		idxinfo=pysam.idxstats(i)
		for j in idxinfo:
			l=(j.strip()).split('\t')
			if l[0]=='*': continue
			if int(l[2])+int(l[3]) > 0: dgdic[l[0]]=i
###########################################################
if not os.path.exists(fastafile):
	usage()
	sys.exit('Fasta file %s not found.' %(fastafile))
if not os.path.exists(fastafile+'.fai'):
	sys.stderr.write('Indexing Fasta file.\n')
	pysam.faidx(fastafile)
###########################################################
# Check reference for name consistency
grefs=dgdic.keys()
rrefs={}
ridxinfo=pysam.idxstats(bamfile)
print ridxinfo
for j in ridxinfo:
	print j
	l=(j.strip()).split('\t')
	if l[0]=='*': continue
	if int(l[2])+int(l[3]) > 0: rrefs[l[0]]=int(l[1])
frefs=[]
fidxinfo=open(fastafile+'.fai')
for j in fidxinfo:
	l=(j.strip()).split('\t')
	if l[0]=='': continue
	frefs.append(l[0])
fidxinfo.close()
#in rna-seq
rnof=[]
for i in rrefs.keys():
	if i not in frefs: sys.stderr.write('WARNING: Region %s in RNA-Seq not found in reference file.\n' %(i))
if len(gbamfile)!=0:
	for i in grefs:
		if i not in frefs: sys.stderr.write('WARNING: Region %s in DNA-Seq not found in reference file.\n' %(i))	
###########################################################

###########################################################
# Annotation file for working regions
if uwf:
	if not os.path.exists(wfile):
		usage()
		sys.exit('GFF file %s not found.' %(wfile))
	if sortann:
		if not whereis('grep'): sys.exit('grep command not found.')
		if not whereis('sort'): sys.exit('sort command not found.')
		sys.stderr.write('Sorting GFF file.\n')
		scmd='grep ^"#" %s; grep -v ^"#" %s | sort -k1,1 -k4,4n > %s' %(wfile,wfile,'workf_%s'%(pid))
		os.system(scmd)
		os.rename(wfile,wfile+'_old')
		os.rename('workf_%s'%(pid),wfile)
	if not os.path.exists(wfile+'.tbi'):
		sys.stderr.write('Indexing GFF file.\n')
		wfile=pysam.tabix_index(wfile, preset='gff')
###########################################################
# Annotation file for strand detection
if uann: 
	getstrand=0
	if not os.path.exists(annfile):
		usage()
		sys.exit('Annotation file %s not found.' %(annfile))
	if sortann:
		if not whereis('grep'): sys.exit('grep command not found.')
		if not whereis('sort'): sys.exit('sort command not found.')
		sys.stderr.write('Sorting annotation file.\n')
		scmd='grep ^"#" %s; grep -v ^"#" %s | sort -k1,1 -k4,4n > %s' %(annfile,annfile,'annotation_%s'%(pid))
		os.system(scmd)
		os.rename(annfile,annfile+'_old')
		os.rename('annotation_%s'%(pid),annfile)
	if not os.path.exists(annfile+'.tbi'):
		sys.stderr.write('Indexing annotation file.\n')
		annfile=pysam.tabix_index(annfile, preset='gff')
###########################################################
# Annotation file to exclude positions
if expos:
	if not os.path.exists(exfile):
		usage()
		sys.exit('File %s not found.' %(exfile))
	if sortann:
		if not whereis('grep'): sys.exit('grep command not found.')
		if not whereis('sort'): sys.exit('sort command not found.')
		sys.stderr.write('Sorting file.\n')
		scmd='grep ^"#" %s; grep -v ^"#" %s | sort -k1,1 -k4,4n > %s' %(exfile,exfile,'exfile_%s'%(pid))
		os.system(scmd)
		os.rename(exfile,exfile+'_old')
		os.rename('exfile_%s'%(pid),exfile)
	if not os.path.exists(exfile+'.tbi'):
		sys.stderr.write('Indexing %s file.\n' %(exfile))
		exfile=pysam.tabix_index(exfile, preset='gff')	
###########################################################
#mainbam=pysam.Samfile(bamfile,"rb")
#regions=mainbam.references
#regionslens=mainbam.lengths
#mainbam.close()
dicregions=dict(rrefs.items())
#dicregions=dict([(regions[x],regionslens[x]) for x in range(len(regions))])
chrs=[x for x in dicregions.keys() if x not in nochrs]
sys.stderr.write('Analysis on %i regions.\n' %(len(chrs)))
###########################################################
if infolder!='': outfolder=os.path.join(outfolder_,'DnaRna_%s_%s' %(infolder,pid))
else: outfolder=os.path.join(outfolder_,'DnaRna_%s' %(pid))
if not os.path.exists(outfolder):
	splitfolder=os.path.split(outfolder)
	if not os.path.exists(splitfolder[0]): os.mkdir(splitfolder[0])
	os.mkdir(outfolder)	
outtable=os.path.join(outfolder,'outTable_%s' %(pid))
if slist:
	slistfile=os.path.join(outfolder,'outPileupRNA_%s' %(pid))
	if len(gbamfile)!=0: gslistfile=os.path.join(outfolder,'outPileupDNA_%s' %(pid))
#write command line and input parameters
f=open(os.path.join(outfolder,'parameters.txt'),'w')
f.writelines(params)
f.close()
###########################################################

def exploreBAM(myinput):
	isgbam=1
	inputs=myinput.split('$')
	chr,bamfile=inputs[0],inputs[1]
	if not dgdic.has_key(chr): isgbam=0
	outfile=os.path.join(outfolder,'table_%s_%s'%(chr,pid))
	if slist:
		if gziptab: outrna=gzip.open(os.path.join(outfolder,'pileupRNA_%s_%s.gz'%(chr,pid)),'wb')
		else: outrna=open(os.path.join(outfolder,'pileupRNA_%s_%s'%(chr,pid)),'w')
		if not nogbam and isgbam:
			if gziptab: outdna=gzip.open(os.path.join(outfolder,'pileupDNA_%s_%s.gz'%(chr,pid)),'wb')
			else: outdna=open(os.path.join(outfolder,'pileupDNA_%s_%s'%(chr,pid)),'w')
	d,di,gd={},{},{}
	bam=pysam.Samfile(bamfile,"rb")
	if not nogbam and isgbam:
		gbam=pysam.Samfile(dgdic[chr],"rb")
	fasta=pysam.Fastafile(fastafile)
	lenregion=dicregions[chr]
	if uann: tabix=pysam.Tabixfile(annfile)
	if expos: extabix=pysam.Tabixfile(exfile)
	if uwf: wtabix=pysam.Tabixfile(wfile)
	if gziptab: out=gzip.open(outfile+'.gz','wb')
	else: out=open(outfile,'w')
	sys.stderr.write('Started analysis on region: %s\n'%(chr))
	if blatr:
		badblat=os.path.join(blatfolder,'blatseqs_%s.bad'%(chr))
		if os.path.exists(badblat):
			sys.stderr.write('Using Blat mapping for region %s\n'%(chr))
			f=open(badblat)
			for i in f:
				l=(i.strip()).split()
				d[l[0]+'_'+l[1]]=int(l[1])
			f.close()
			sys.stderr.write('Found %i reads for region %s\n'%(len(d),chr))
	if gblatr:
		gbadblat=os.path.join(gblatfolder,'blatseqs_%s.bad'%(chr))
		if os.path.exists(gbadblat):
			sys.stderr.write('Using Blat mapping for DNA region %s\n'%(chr))
			f=open(gbadblat)
			for i in f:
				l=(i.strip()).split()
				gd[l[0]+'_'+l[1]]=int(l[1])
			f.close()
			sys.stderr.write('Found %i reads for region %s\n'%(len(gd),chr))
	if exss:
		if os.path.exists(splicefile):
			sys.stderr.write('Loading known splice sites for region %s\n'%(chr))
			f=open(splicefile)
			for i in f:
				l=(i.strip()).split()
				if l[0]!=chr: continue
				st,tp,cc=l[4],l[3],int(l[1])
				if st=='+' and tp=='D':
					for j in range(nss): di[cc+(j+1)]=0
				if st=='+' and tp=='A':
					for j in range(nss): di[cc-(j+1)]=0
				if st=='-' and tp=='D': 	
					for j in range(nss): di[cc-(j+1)]=0
				if st=='-' and tp=='A':
					for j in range(nss): di[cc+(j+1)]=0	
			f.close()
			sys.stderr.write('Loaded %i positions for %s\n'%(len(di),chr))
	for kpos in range(0,lenregion,chunckval):
		startk,endk=kpos,(kpos+chunckval)-1
		# check features in the give region id GFF provided
		if uwf:
			if chr in wtabix.contigs:
				wfeat=[parseFeat(feat) for feat in wtabix.fetch(reference=chr,start=startk,end=endk)]
				if len(wfeat)==0: continue
				wcoords,startk,endk=newCoords(wfeat,startk,endk) 
			else: continue
		# explore dna-seq bam
		#####################
		gdic={}
		if not nogbam and isgbam:
			for pileupcolumn in gbam.pileup(chr,startk,endk):
				if uwf and not checkPos(wcoords,pileupcolumn.pos): continue
				if not startk<=pileupcolumn.pos<=endk: continue
				gref=fasta.fetch(chr,pileupcolumn.pos,pileupcolumn.pos+1).upper()
				gseq,gqual,gstrand,gsqual,gblatc='',0,'','',''
				if grmsh:
					if ((pileupcolumn.pos+1)-ghomo)-1 < 0: sequp=''
					else: sequp=(fasta.fetch(chr,((pileupcolumn.pos+1)-ghomo)-1,(pileupcolumn.pos+1)-1)).upper()
					seqdw=(fasta.fetch(chr,pileupcolumn.pos+1,(pileupcolumn.pos+1)+ghomo)).upper()
				for pileupread in pileupcolumn.pileups: # per ogni base dell'allineamento multiplo
					gs,gq,gt,gqq=pileupread.alignment.seq[pileupread.qpos].upper(),ord(pileupread.alignment.qual[pileupread.qpos])-gQVAL,'*',pileupread.alignment.qual[pileupread.qpos]
					# multiple hits
					if gexh and pileupread.alignment.is_secondary: continue
					# duplicates
					if gexd and pileupread.alignment.is_duplicate: continue
					# se paired end
					if gconc and pileupread.alignment.is_paired:
						# se non concordanti
						if not pileupread.alignment.is_proper_pair: continue
						# se concordanti ma nello stesso orientamento
						flag=pileupread.alignment.flag
						if pileupread.alignment.is_duplicate: flag=flag-1024
						if pileupread.alignment.is_secondary: flag=flag-256
						if flag in [67,131,115,179]: continue
					# mapping quality
					if gmq and pileupread.alignment.mapq < gMAPQ: continue  
					#se la qualita' >= della qualita' minima
					if gq >= gMQUAL and pileupcolumn.pos in pileupread.alignment.positions:	
						if grmnuc:
							#grlen=pileupread.alignment.qlen #lunghezza della specifica read
							#print rlen,pileupread.qpos,pileupread.alignment.qstart,pileupread.alignment.qend
							# verifica se il nuc deve essere rimosso alle estremita' nel range x-y
							# testare il forward
							#gqp=pileupread.qpos-pileupread.alignment.qstart
							#print qp,pileupread.alignment.is_reverse
							#print (rlen-qp)-1
							#if pileupread.alignment.is_reverse:
							#	if (grlen-gqp)-1 < grmp[0]:continue
							#	if (grlen-gqp)-1 > ((grlen)-grmp[1])-1: continue
							#else:
							#	if gqp<grmp[0]:continue
							#	if gqp>(grlen-grmp[1])-1: continue
							grlen=pileupread.alignment.rlen #pileupread.alignment.qlen #lunghezza della specifica read
							gqp=pileupread.qpos #pileupread.qpos-pileupread.alignment.qstart
							if pileupread.alignment.is_reverse:
								if gqp>(grlen-grmp[0])-1: continue
								if gqp<grmp[1]:continue
							else:
								if gqp<grmp[0]:continue
								if gqp>(grlen-grmp[1])-1: continue
						# se la read di appartenenza non mappa in modo univoco con Blat
						if gblatr:
							rt=0
							if pileupread.alignment.is_read1: rt=1
							elif pileupread.alignment.is_read2: rt=2
							rname=pileupread.alignment.qname+'_%i'%(rt)
							if gd.has_key(rname): gblatc+='0' #continue
							else: gblatc+='1'
						# se la base e' diversa dal reference
						# se in regione omopolimerica scarta
						if grmsh and rmHomo(sequp,seqdw,ghomo,gref): continue
						gseq+=gs
						gqual+=gq
						gstrand+=gt
						gsqual+=gqq
				if gseq.strip()!='':
					if gblatr:
						if testBlat(gblatc): gseq,gqual,gsqual,gstrand=normByBlat(gseq,gstrand,gsqual,gblatc,gQVAL)
						else: continue					
					gcov,gbcomp,gsubs,gfreq=BaseCount(gseq,gref,gmmf,0)
					if gcov < gMINCOV: continue
					gmqua=meanq(gqual,len(gseq))
					ghinfo=0 # non omozigote
					if gsubs=='-': ghinfo=1 # omozigote
					gdic[pileupcolumn.pos]=([str(gcov),gmqua,str(gbcomp),gsubs,gfreq],ghinfo)
					if slist:
						if not nogbam and isgbam: outdna.write('\t'.join([chr,str(pileupcolumn.pos+1),gref,gseq,gsqual])+'\n')
		#####################
		# explore rna-seq bam
		for pileupcolumn in bam.pileup(chr,startk,endk):
			if uwf and not checkPos(wcoords,pileupcolumn.pos): continue
			if not startk<=pileupcolumn.pos<=endk: continue
			ref=fasta.fetch(chr,pileupcolumn.pos,pileupcolumn.pos+1).upper()
			seq,qual,strand,squal,blatc='',0,'','',''
			if rmsh:
				if ((pileupcolumn.pos+1)-homo)-1 < 0: sequp=''
				else: sequp=(fasta.fetch(chr,((pileupcolumn.pos+1)-homo)-1,(pileupcolumn.pos+1)-1)).upper()
				seqdw=(fasta.fetch(chr,pileupcolumn.pos+1,(pileupcolumn.pos+1)+homo)).upper()
			for pileupread in pileupcolumn.pileups: # per ogni base dell'allineamento multiplo
				s,q,t,qq=pileupread.alignment.seq[pileupread.qpos].upper(),ord(pileupread.alignment.qual[pileupread.qpos])-QVAL,'*',pileupread.alignment.qual[pileupread.qpos]
				# escludi posizioni introniche nei pressi di splice sites
				if exss and di.has_key(pileupcolumn.pos+1): continue
				# multiple hit
				if exh and pileupread.alignment.is_secondary: continue
				# duplicates
				if exd and pileupread.alignment.is_duplicate: continue
				# se paired end
				if conc and pileupread.alignment.is_paired:
					# se non concordanti
					if not pileupread.alignment.is_proper_pair: continue
					# se concordanti ma nello stesso orientamento
					flag=pileupread.alignment.flag
					if pileupread.alignment.is_duplicate: flag=flag-1024
					if pileupread.alignment.is_secondary: flag=flag-256
					if flag in [67,131,115,179]: continue
				# mapping quality
				if mq and pileupread.alignment.mapq < MAPQ: continue  
				#se la qualita' >= alla qualita' minima
				if q >= MQUAL and pileupcolumn.pos in pileupread.alignment.positions:
					#tags=dict(pileupread.alignment.tags)
					#deduci la strand per ogni posizione
					if getstrand:
						#usa le info del mapping se strand oriented 
						if pileupread.alignment.is_read1:
							if unchange1:
								if pileupread.alignment.is_reverse: t='-'
								else: t='+'
							else:
								if pileupread.alignment.is_reverse: t='+'
								else: t='-'
						elif pileupread.alignment.is_read2:
							if unchange2:
								if pileupread.alignment.is_reverse: t='-'
								else: t='+'
							else:
								if pileupread.alignment.is_reverse: t='+'
								else: t='-'
						else: # for single ends
							if unchange1:
								if pileupread.alignment.is_reverse: t='-'
								else: t='+'
							else:
								if pileupread.alignment.is_reverse: t='+'
								else: t='-'																					
					if rmnuc:
						#rlen=pileupread.alignment.qlen #lunghezza della specifica read
						#print rlen,pileupread.qpos,pileupread.alignment.qstart,pileupread.alignment.qend
						# verifica se il nuc deve essere rimosso alle estremita' nel range x-y
						# testare il forward
						#qp=pileupread.qpos-pileupread.alignment.qstart
						#print pileupread.qpos,pileupread.alignment.rlen
						#print (rlen-qp)-1
						#if pileupread.alignment.is_reverse:
						#	if (rlen-qp)-1 < rmp[0]:continue
						#	if (rlen-qp)-1 > ((rlen)-rmp[1])-1: continue
						#else:
						#	if qp<rmp[0]:continue
						#	if qp>(rlen-rmp[1])-1: continue				
						rlen=pileupread.alignment.rlen #pileupread.alignment.qlen #lunghezza della specifica read
						qp=pileupread.qpos #pileupread.qpos-pileupread.alignment.qstart
						if pileupread.alignment.is_reverse:
							if qp>(rlen-rmp[0])-1: continue
							if qp<rmp[1]:continue
						else:
							if qp<rmp[0]:continue
							if qp>(rlen-rmp[1])-1: continue
					# se la read di appartenenza non mappa in modo univoco con Blat
					if blatr:
						rt=0
						if pileupread.alignment.is_read1: rt=1
						elif pileupread.alignment.is_read2: rt=2
						rname=pileupread.alignment.qname+'_%i'%(rt)
						if d.has_key(rname): blatc+='0' #continue
						else: blatc+='1'
					# se la base e' diversa dal reference
					# se in regione omopolimerica scarta
					if rmsh and rmHomo(sequp,seqdw,homo,ref): continue
					seq+=s
					qual+=q
					strand+=t
					squal+=qq
			if seq.strip()!='':
				if blatr:
					if testBlat(blatc): seq,qual,squal,strand=normByBlat(seq,strand,squal,blatc,QVAL)
					else: continue
				mystrand='2'
				#print seq,strand,strand.count('+'),strand.count('-')
				if uann and not getstrand:
					if chr in tabix.contigs:
						sres=[kk.strand for kk in tabix.fetch(reference=chr,start=(pileupcolumn.pos),end=(pileupcolumn.pos+1),parser=pysam.asGTF())]
						mystrand=vstrand(sres)
				if getstrand and not uann:
					mystr=vstand(strand)
					if mystr=='-': mystrand='0'
					elif mystr=='+': mystrand='1'
					else: mystrand='2'
				if mystrand=='0':
					seq=comp(seq)
					ref=comp(ref)
				if mystrand in ['0','1'] and corrstr:
					seq,qual,squal=normByStrand(seq,strand,squal,mystrand)
				cov,bcomp,subs,freq=BaseCount(seq,ref,mmf,vnuc)
				if cov < MINCOV: continue
				if exms and subs.count(' ')>0: continue
				mqua=meanq(qual,len(seq))
				if expos:
					if chr in extabix.contigs:
						exres=[kk for kk in extabix.fetch(reference=chr,start=(pileupcolumn.pos),end=(pileupcolumn.pos+1))]
						if len(exres)>0: continue
				# se la sostituzione non e' in usubs
				if exinv and subs=='-': continue
				if not checkSubs(subs): continue
				#print out rna-seq info + dna-seq
				if gdic.has_key(pileupcolumn.pos): # abbiamo l'informazione genomica
					if exnonh and not gdic[pileupcolumn.pos][1]: continue
					if mystrand=='0':
						gdic[pileupcolumn.pos][0][2]=comp2(eval(gdic[pileupcolumn.pos][0][2]))
						gdic[pileupcolumn.pos][0][3]=comp(gdic[pileupcolumn.pos][0][3])
					line='\t'.join([chr,str(pileupcolumn.pos+1),ref,mystrand,str(cov),mqua,str(bcomp),subs,freq]+gdic[pileupcolumn.pos][0])+'\n'
					out.write(line)
				else:
					if exnosupp: continue
					line='\t'.join([chr,str(pileupcolumn.pos+1),ref,mystrand,str(cov),mqua,str(bcomp),subs,freq]+['-','-','-','-','-'])+'\n'
					out.write(line)
				if slist: outrna.write('\t'.join([chr,str(pileupcolumn.pos+1),ref,seq,squal])+'\n')
	bam.close()
	if not nogbam and isgbam: gbam.close()
	fasta.close()
	out.close()
	if uwf: wtabix.close()
	if uann: tabix.close()
	if expos: extabix.close()
	if slist:
		outrna.close()
		if not nogbam and isgbam: outdna.close()
	sys.stderr.write('Job completed for region: %s\n'%(chr))
	
def do_work(q):
	while True:
		try:
			x=q.get(block=False)
			exploreBAM(x)
		except Empty:
			break

work_queue = Queue()
for i in chrs:
	strinput=i+'$'+bamfile
	work_queue.put(strinput)
processes=[Process(target=do_work, args=(work_queue,)) for i in range(NCPU)]
for t in processes:
	t.start()
for t in processes:
	t.join()
time.sleep(0.5)
#
head='Region\tPosition\tReference\tStrand\tCoverage-q%i\tMeanQ\tBaseCount[A,C,G,T]\tAllSubs\tFrequency\tgCoverage-q%i\tgMeanQ\tgBaseCount[A,C,G,T]\tgAllSubs\tgFrequency\n' %(MQUAL,gMQUAL)
sys.stderr.write('Merging Tables.\n')
if gziptab: o=gzip.open(outtable+'.gz','wb')
else: o=open(outtable,'w')
o.write(head)
if slist:
	if gziptab: o2=gzip.open(slistfile+'.gz','wb')
	else: o2=open(slistfile,'w')
	if len(gbamfile)!=0:
		if gziptab: o3=gzip.open(gslistfile+'.gz','wb')
		else: o3=open(gslistfile,'w')
for i in chrs:
	if gziptab: tabfile=os.path.join(outfolder,'table_%s_%s.gz' %(i,pid))
	else: tabfile=os.path.join(outfolder,'table_%s_%s' %(i,pid))
	if os.path.exists(tabfile):
		if gziptab: f=gzip.open(tabfile,'rb')
		else: f=open(tabfile)
		for j in f: o.write(j)
		f.close()
		os.remove(tabfile)
	if slist:
		if len(gbamfile)!=0:
			if gziptab: dnafile=os.path.join(outfolder,'pileupDNA_%s_%s.gz' %(i,pid))
			else: dnafile=os.path.join(outfolder,'pileupDNA_%s_%s' %(i,pid))
			if os.path.exists(dnafile):
				if gziptab: f=gzip.open(dnafile,'rb')
				else: f=open(dnafile)
				for j in f: o3.write(j)
				f.close()
				os.remove(dnafile)
		if gziptab: rnafile=os.path.join(outfolder,'pileupRNA_%s_%s.gz' %(i,pid))
		else: rnafile=os.path.join(outfolder,'pileupRNA_%s_%s' %(i,pid))
		if os.path.exists(rnafile):
			if gziptab: f=gzip.open(rnafile,'rb')
			else: f=open(rnafile)
			for j in f: o2.write(j)
			f.close()
			os.remove(rnafile) 
o.close()
if slist:
	o2.close()
	if len(gbamfile)!=0: o3.close()
sys.stderr.write('Results saved on %s\n'%(outtable))
if slist:
	if len(gbamfile)!=0: sys.stderr.write('Pileup for DNA saved on %s\n'%(gslistfile))
	sys.stderr.write('Pileup for RNA saved on %s\n'%(slistfile))
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stderr.write("Script time --> END: %s\n"%(script_time))

