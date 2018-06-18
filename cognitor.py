#!/usr/bin/python

import sys
sys.dont_write_bytecode = True

import os, shutil
import subprocess 
from query2subjectFix import createFile
from config import cognitor_settings

try:

	## deletes old files then recreates them ##
	folders = [cognitor_settings['main_directory'] + 'BLASTss', cognitor_settings['main_directory'] + 'BLASTno', cognitor_settings['main_directory'] + 'BLASTff', cognitor_settings['main_directory'] + 'BLASTcogn', cognitor_settings['main_directory'] + 'BLASTconv', cognitor_settings['main_directory']];
	for x in folders:

		try:
			os.remove(x)
		except OSError as e: 
		    pass;

	if os.path.exists(cognitor_settings['main_directory']):
		shutil.rmtree(cognitor_settings['main_directory'])

	if not os.path.exists(cognitor_settings['main_directory']):
	    os.makedirs(cognitor_settings['main_directory'])
	    os.makedirs(cognitor_settings['main_directory'] + "BLASTss")
	    os.makedirs(cognitor_settings['main_directory'] + "BLASTno")
	    os.makedirs(cognitor_settings['main_directory'] + "BLASTff")
	    os.makedirs(cognitor_settings['main_directory'] + "BLASTcogn")
	    os.makedirs(cognitor_settings['main_directory'] + "BLASTconv")

	createFile(cognitor_settings['query_csv_file'], cognitor_settings['subject_fa']);

	print('files reset');

	#makes databases 
	os.system("makeblastdb -logfile " + cognitor_settings['main_directory'] + "makeblastdb_sub.log -in " + cognitor_settings['subject_fa'] +" -dbtype prot -out " + cognitor_settings['query_db'] )
	os.system("makeblastdb -logfile " + cognitor_settings['main_directory'] + "makeblastdb_cog.log -in " + cognitor_settings['cog_domains_fa'] + "  -dbtype prot -out " + cognitor_settings['cogs_db'] )
	print('Databases Created');

	#creates file that bridges together query and subject sequences, then it hashes the file 
	os.system("cat " + cognitor_settings['query_csv_file'] + " " + cognitor_settings['cog_domains'] + " > "+ cognitor_settings['query_plus_subject_csv'] )
	os.system("COG_software/COGmakehash/COGmakehash -i="+ cognitor_settings['query_plus_subject_csv'] +" -o=" + cognitor_settings['main_directory'] + "BLASTcogn -s='' -n=1") 
	print('Hash Files created')

	os.system("psiblast -query " + cognitor_settings['subject_fa'] + " -db " + cognitor_settings['query_db'] + " -show_gis -outfmt 7 -num_threads " + str(cognitor_settings['num_threads']) + " -num_alignments 10 -dbsize 100000000 -comp_based_stats F -seg no -out " + cognitor_settings['main_directory'] + "BLASTss/QuerySelf.tab")
	print('blast self vs self finished');
	os.system("psiblast -query " + cognitor_settings['subject_fa'] + " -db " + cognitor_settings['cogs_db'] + " -show_gis -outfmt 7 -num_threads " + str(cognitor_settings['num_threads']) + " -num_alignments 1000 -dbsize 100000000 -comp_based_stats F -seg no -out " + cognitor_settings['main_directory'] + "BLASTno/QueryCOGs.tab") 
	print('blast all vs all unfiltered finished');
	os.system("psiblast -query " + cognitor_settings['subject_fa'] + " -db " + cognitor_settings['cogs_db'] + " -show_gis -outfmt 7 -num_threads " + str(cognitor_settings['num_threads']) + " -num_alignments 1000 -dbsize 100000000 -comp_based_stats T -seg yes -out " + cognitor_settings['main_directory'] + "BLASTff/QueryCOGs.tab")
	print('blast all vs all filtered finished')

	print("cognitor started");

	os.system("COG_software/COGreadblast/COGreadblast -d=" + cognitor_settings['main_directory'] + "BLASTcogn -u=" + cognitor_settings['main_directory'] + "BLASTno -f=" + cognitor_settings['main_directory'] + "BLASTff -s=" + cognitor_settings['main_directory'] + "BLASTss -e=0.1 -q=2 -t=2")
	os.system("COG_software/COGcognitor/COGcognitor -i=" + cognitor_settings['main_directory'] + "BLASTcogn -t="+cognitor_settings['cog_domains']+" -q=" + cognitor_settings['query_csv_file'] + " -o=" + cognitor_settings['result_file'] )

	shutil.move("cognitor.log", cognitor_settings['main_directory'] + "cognitor.log")
	shutil.move("conflict.txt", cognitor_settings['main_directory'] + "conflict.txt")

	print("cognitor finished");

	pass

except Exception as e:

	print('Error: ', str(e))

