
main_directory = '../single_gene_results/gene_set_test_3/data/';
result_directory = '../single_gene_results/gene_set_test_3/';
genomes_directory = '../genomes/';


cognitor_settings = {
	
	#main main_directory for output and results 
	"main_directory" : main_directory,
	"genomes_directory" : genomes_directory,

	#hash and combined query csv temp files 
	#These are required by the cog software 
	"query_csv_file" : main_directory + "query.p2o.csv",
	"query_plus_subject_csv" : main_directory + "tmp.p2o.csv",

	#query set
	"subject_fa" : genomes_directory + "gene_set.fa",

	#cogs to test against 
	#ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cog2003-2014.csv
	"cog_domains" : genomes_directory + "cogs.p2o.csv",
	# ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/prot2003-2014.fa.gz
	"cog_domains_fa" : genomes_directory + "cogs.fa",

	#Blast Databases
	"cogs_db" : main_directory + "COGs",
	"query_db" : main_directory + "GenThree",

	#system options 
	"num_threads" : 8,

	#the result of your analysis 
	"result_file" : result_directory + "results.COG.csv"
	
}

main_directory = '../test2/triTemp/';

cog_maker_settings = {
	
	#main main_directory for output and results 
	"main_directory" : main_directory,
	"genomes_directory" : genomes_directory, 

	#hash and combined query csv temp files 
	#These are required by the cog software 
	"query_csv_file" : main_directory + "query.p2o.csv",

	#blast Database
	"query_db" : main_directory + "GenThree",

	#query set
	"subject_fa" : genomes_directory + "GenThree.fa",

	# job file 
	"job_file" : genomes_directory + "job_file.csv",

	#lineage specific expansion file
	"lse_output_file" : main_directory + "GenThree.cls.csv",

	#final results file 
	"result_file" : main_directory + "results.cls.csv",

	#system options 
	"num_threads" : 8,
	
}
