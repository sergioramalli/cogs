#!/usr/bin/python

import sys
sys.dont_write_bytecode = True

import os, shutil
import subprocess 
from query2subjectFix import createFile
from config import cog_maker_settings


# $ cat Gen1.fa Gen2.fa Gen3.fa > GenThree.fa # pools the sets together

try:

	## deletes old files then recreates them ##
	try:
		os.remove(cog_maker_settings['main_directory'])
	except OSError as e: 
	    pass;

	if os.path.exists(cog_maker_settings['main_directory']):
		shutil.rmtree(cog_maker_settings['main_directory'])
		
	if os.path.exists("./data/tmp/"):
		shutil.rmtree("./data/tmp/")

	if not os.path.exists(cog_maker_settings['main_directory']):
	    os.makedirs("./data/tmp/")
	    os.makedirs(cog_maker_settings['main_directory'])
	    os.makedirs(cog_maker_settings['main_directory'] + "BLASTno")
	    os.makedirs(cog_maker_settings['main_directory'] + "BLASTff")
	    os.makedirs(cog_maker_settings['main_directory'] + "BLASTconv")
	    # BLASTconv

	createFile(cog_maker_settings['query_csv_file'], cog_maker_settings['subject_fa']);
	print('files reset');

	# formats a BLASTable database
	os.system("makeblastdb -in " + cog_maker_settings['subject_fa'] + " -dbtype prot -out " + cog_maker_settings['query_db'])

	# unfiltered blast results in the blastno directory
	os.system("psiblast -query " + cog_maker_settings['subject_fa'] + " -db " + cog_maker_settings['query_db'] + " -show_gis -outfmt 7 -num_threads " + str(cog_maker_settings['num_threads']) + " -num_alignments 1000 -dbsize 100000000 -comp_based_stats F -seg no -out "+ cog_maker_settings['main_directory'] +"BLASTno/ThreeByThree.tab ")

	# $ # filtered BLAST results in the ./BLASTff/ directory
	os.system("psiblast -query " + cog_maker_settings['subject_fa'] + " -db " + cog_maker_settings['query_db'] + " -show_gis -outfmt 7 -num_threads " + str(cog_maker_settings['num_threads']) + " -num_alignments 1000 -dbsize 100000000 -comp_based_stats T -seg yes -out "+ cog_maker_settings['main_directory'] +"BLASTff/ThreeByThree.tab ")

	os.system("COG_software/COGmakehash/COGmakehash -i=" + cog_maker_settings['query_csv_file'] + " -o="+ cog_maker_settings['main_directory'] +"BLASTconv -s=',' -n=1")

	os.system("COG_software/COGreadblast/COGreadblast -d=" + cog_maker_settings['main_directory'] + "BLASTconv -u=" + cog_maker_settings['main_directory'] + "BLASTno -f=" + cog_maker_settings['main_directory'] + "BLASTff -s=" + cog_maker_settings['main_directory'] + "BLASTno -e=0.1 -q=2 -t=2");

	# create gen file
	os.system("COG_software/COGlse/COGlse -d=" + cog_maker_settings['main_directory'] + "BLASTconv -j=" + cog_maker_settings['job_file'] + " -p=" + cog_maker_settings['query_csv_file'] + " -o=" + cog_maker_settings['lse_output_file'] )

	os.system("COG_software/COGtriangles/COGtriangles -i=" + cog_maker_settings['main_directory'] + "BLASTconv -q=" + cog_maker_settings['query_csv_file'] + " -l=" + cog_maker_settings['lse_output_file'] + " -o="+ cog_maker_settings['result_file'] +" -t=0.5 -e=0.01 -n='CLS' -s=1")

	pass

except Exception as e:

	exc_type, exc_obj, exc_tb = sys.exc_info()
	fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
	print(exc_type, fname, exc_tb.tb_lineno)

