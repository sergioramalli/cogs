#!/usr/bin/python
import os
import subprocess 
from query2subjectFix import createFile

directory = '../cogTemp/';

## Try to delete the file ##
folders = [directory + 'BLASTss', directory + 'BLASTno', directory + 'BLASTff', directory + 'BLASTcogn', directory + 'BLASTconv', directory];
for x in folders:

	try:
		os.remove(x)
	except OSError, e:  ## if failed, report it back to the user ##
	    pass;


if not os.path.exists(directory):
    os.makedirs(directory)

print('files reset');

outputFile = directory + "query.p2o.csv";
inputFile = "../genomes/reduced.fa";

cog_fa = "../genomes/cogs.fa"
cogs_db1 = directory + "COGs"

genThree_fa = "../genomes/reduced.fa";
genThree_db = directory + "GenThree";

cog_domains = '../genomes/cogs.p2o.csv';

createFile(outputFile, inputFile);

# makeblastdb -in GenThree.fa -dbtype prot -out GenThree

os.system("makeblastdb -in " + genThree_fa +" -dbtype prot -out " + genThree_db )
# print(os.system("makeblastdb -in " + cog_fa + "  -dbtype prot -out " + cogs_db1 ))

print('Databases Created');

os.system("cat " + outputFile + " " + cog_domains + " > tmp.p2o.csv" )



#  -in GenThree.fa -dbtype prot -out GenThree

# makeblastdb -in cogs.fa -dbtype prot -out COGs 

# psiblast -query GenThree.fa -db GenThree -show_gis -outfmt 8 -num_threads 4 -num_alignments 10 -dbsize 100000000 -comp_based_stats F -seg no -out BLASTss/QuerySelf.tab 

# psiblast -query GenThree.fa -db COGs -show_gis -outfmt 8 -num_threads 4 -num_alignments 1000 -dbsize 100000000 -comp_based_stats F -seg no -out BLASTno/QueryCOGs.tab 

# psiblast -query GenThree.fa -db COGs -show_gis -outfmt 8 -num_threads 4 -num_alignments 1000 -dbsize 100000000 -comp_based_stats T -seg yes -out BLASTff/QueryCOGs.tab 

# date
# COG_software/COGmakehash/COGmakehash -i=tmp.p2o.csv -o=./BLASTcogn -s="," -n=1 

# COG_software/COGreadblast/COGreadblast -d=./BLASTcogn -u=./BLASTno -f=./BLASTff -s=./BLASTss -e=0.1 -q=2 -t=2

# COG_software/COGcognitor/COGcognitor -i=./BLASTcogn -t=cogs.p2o.csv -q=query.p2o.csv -o=GenQuery.COG.csv 

# date

# ls




# print( subprocess.check_output(['ls', '-l']) ); 