#!usr/bin

#Usage:Correct_strand.py -a /xtdisk/apod-gen/zhutt/reference/Homo_sapiens.sorted.GRCh38.99.gtf.gz -i $1/parallel_table.txt -o $1/parallel_table_corrected.txt

import sys, os, getopt, random, time
try: import pysam
except: sys.exit('Pysam module not found.')
pid=str(os.getpid()+random.randint(0,999999999))

def usage():
	print """
USAGE: python AnnotateTable.py [options]
Options:
-a		Sorted Annotation file
-i		Annotate a file of positions [column1=region, column2=coordinate (1 based)]
		or a single position [region:coordinate (1 based)]
-o		Save lines to a file
-h		Print this help
"""

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:a:o:h',["help"])
except getopt.GetoptError, err:
	print str(err) 
	usage()
	sys.exit()

if len(opts)==0:
	usage()
	sys.exit()

tablefile,outfile,annfile='','',''
for o,a in opts:
	if o in ("-h","--help"):
		usage()
		sys.exit()
	elif o == "-i":
		tablefile = a
	elif o == "-o":
		outfile = a
	elif o == "-a":
		annfile = a
		if annfile=='':
			usage()
			sys.exit('Sorted annotation file not found.')
	else:
		assert False, "unhandled option"

def parse(res):
	d={'+':{},'-':{}}
	anns='+'
	for i in res:
		if i[3]=='+':
			if d['+'].has_key(i[1]):
				if i[0] not in d['+'][i[1]][0]: d['+'][i[1]][0]=d['+'][i[1]][0]+','+i[0]
				if i[2]+'-'+i[0] not in d['+'][i[1]][1]: d['+'][i[1]][1]=d['+'][i[1]][1]+','+i[2]+'-'+i[0]
			else:
				d['+'][i[1]]=[i[0],i[2]+'-'+i[0]]
		elif i[3]=='-':
			if d['-'].has_key(i[1]):
				if i[0] not in d['-'][i[1]][0]: d['-'][i[1]][0]=d['-'][i[1]][0]+','+i[0]
				if i[2]+'-'+i[0] not in d['-'][i[1]][1]: d['-'][i[1]][1]=d['-'][i[1]][1]+','+i[2]+'-'+i[0]
			else:
				d['-'][i[1]]=[i[0],i[2]+'-'+i[0]]
	gip='$'.join(d['+'].keys())
	featp='$'.join([d['+'][x][0] for x in d['+'].keys()])
	tip='$'.join([d['+'][x][1] for x in d['+'].keys()])
	gim='$'.join(d['-'].keys())
	featm='$'.join([d['-'][x][0] for x in d['-'].keys()])
	tim='$'.join([d['-'][x][1] for x in d['-'].keys()])
	p=[featp,gip,tip]
	m=[featm,gim,tim]
	pm=[(featp+'&'+featm).strip('&'),(gip+'&'+gim).strip('&'),(tip+'&'+tim).strip('&')]
	if len(d['+'])==0 and len(d['-'])!=0: anns='-'
	if len(d['+'])==0: p=['-','-','-']
	if len(d['-'])==0: m=['-','-','-']
	if len(d['+'])==0 and len(d['-'])==0:
		pm=['-','-','-']
		anns='+-'
	if len(d['+'])!=0 and len(d['-'])!=0: anns='+-'
	return (p,m,pm,anns)

script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stderr.write("Script time --> START: %s\n" %(script_time))


if not os.path.exists(annfile+'.tbi'):
	sys.stderr.write('Indexing %s file.\n' %(annfile))
	annfile=pysam.tabix_index(annfile, preset='gff')

tabix=pysam.Tabixfile(annfile)
contig=tabix.contigs
o=open(outfile,'w')
f=open(tablefile)
for i in f:
	if i.strip()=='': continue
	if i.startswith('Region'):
		h=[i.strip()]
		o.write('\t'.join(h)+'\n')
	if i.startswith('#'): continue
	l=(i.strip()).split('\t')
	chr,pos=l[0],int(l[1])-1
	sres=[]
	if chr in contig:
		sres=[(kk.feature,kk.gene_id,kk.gene_name,kk.strand) for kk in tabix.fetch(reference=chr,start=pos,end=pos+1,parser=pysam.asGTF())]
		ann=parse(sres) #(p,m,pm,anns)
		if ann[3]=='+-': l[3]=2
		elif ann[3]=='+': l[3]=1
		elif ann[3]=='-':l[3]=0
		o.write('\t'.join('%s' %id for id in l)+'\n')
tabix.close()
o.close()
sys.stderr.write("Table saved on %s\n" %(outfile))
script_time=time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
sys.stderr.write("Script time --> END: %s\n" %(script_time))
