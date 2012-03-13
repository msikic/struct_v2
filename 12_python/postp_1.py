#!/usr/bin/env python
import sys, os, getopt, subprocess
from multiprocessing import Pool

def usage():
	print "Usage: postp_1.py -s struct_out -p parameters_file -t number_of_used_lines"

def main():                                        

# parsing input values
    
	top_number = ""
	db_search = ""
	params_file = ""
		
	try:                                
		opts, args = getopt.getopt(sys.argv[1:], "hs:p:t:", ["help", "struct=", "param=", "top="]) 
	except getopt.getopt.error, msg:           
		print msg                          
		sys.exit(2) 
	for opt, arg in opts:
		if opt in ("-h", "--help"):      
			usage()                     
			sys.exit() 
		elif opt in ("-s", "--struct"):      
			db_search = arg
		elif opt in ("-p", "--param"):      
			params_file = arg
		elif opt in ("-t", "--top"):      
			top_number = arg
	
# options - should be changed in future 	
	struct    = "/home/miles/Projects/03_struct/struct";
	dbdir     = "dbis"; 
	pdbdir    = "pdbs";
	pdbaffine = "/home/miles/Projects/03_struct/15_postproc/perlscr/pdb_affine_tfm.pl";
	pdb_extr_chain = "/home/miles/Projects/03_struct/15_postproc/perlscr/pdb_extract_chain.pl";

# checking if all directories, scripts and files exist
	necessary_data = [struct, dbdir, pdbdir, pdbaffine, pdb_extr_chain]
	for data in necessary_data:
		if not os.path.exists(data):
			print "{0} path doesn't exist".format(data)
			sys.exit()

	directories = ["tfms", "maps"]
	for d in directories:
		if not os.path:
			os.mkdir(d)
# calling sort and sorting db_search file using the 5th column

	p1 = subprocess.Popen("sort -grk 5 {0}".format(db_search), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	p2 = subprocess.Popen("head -n {0}".format(top_number), shell=True, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	top_list = p2.stdout.read().split('\n')
	top_list=filter(lambda x: len(x)>0, top_list) 

	for hit in top_list:
		hit_list = hit.split(None, 2)
		
		qry_name = hit_list[0]
		qry_pdb_name = qry_name[0:4]		
		qry_db = "".join([dbdir, "/", qry_name, ".db"])	
		qry_pdb = "".join([pdbdir, "/", qry_pdb_name, ".pdb"])		
		qry_chain = qry_name[4:5]


		tgt_name = hit_list[1]
		tgt_pdb_name = tgt_name[0:4]
		tgt_db = "".join([dbdir, "/", tgt_name, ".db"])	
		tgt_pdb = "".join([pdbdir, "/", tgt_pdb_name, ".pdb"])	
		tgt_chain = tgt_name[4:5]

#print tgt_db

		for file_name in [tgt_db, tgt_pdb, qry_db, qry_pdb]:
			if not os.path.exists(file_name):
				print "{0} file doesn't exist".format(file_name)
				sys.exit()
			if not os.path.getsize(file_name):
				print "{0} file is empty".format(file_name)
				sys.exit()
			
# creating a new paramters file and adding postprecesing info

		params_file_new_name = "".join([params_file, tgt_name])
		cmd = "grep -v pdf_ {0} | grep -v chain_ | grep -v postp | grep -v verbose > {1}".format(params_file, params_file_new_name)
		os.system(cmd)

		params_file_new = open(params_file_new_name, 'a')
		params_file_new.write("postp\n")
		params_file_new.write("".join(["pdbf_tgt ", tgt_pdb, "\n"]))
		params_file_new.write("".join(["chain_tgt ", tgt_chain, "\n"]))
		params_file_new.write("".join(["pdbf_qry ", qry_pdb, "\n"]))
		params_file_new.write("".join(["chain_qry ", qry_chain, "\n"]))
		params_file_new.write("verbose\n")
		params_file_new.close()

# run struct

		outname = "_".join([qry_name, tgt_name])
		struct_out = "".join([outname, ".long_out"])
		struct_tfm = "".join([outname, ".struct_out.for_postp"])
		cmd = " ".join([struct, tgt_db, qry_db, params_file_new_name, ">", struct_out])
		os.system(cmd)
	
#os.system("grep -v done {0}.struct_out >> postp_out".format(outname))
		os.rename("{0}.struct_out.for_postp".format(outname), "tfms/{0}_to_{1}.tfm".format(tgt_name, qry_name))
		os.rename(struct_out, "maps/{0}_to_{1}.map".format(tgt_name, qry_name))

		os.remove(params_file_new_name)

	for hit in top_list:
		hit_list = hit.split(None, 2)
		
		qry_name = hit_list[0]
		tgt_name = hit_list[1]
		outname = "_".join([qry_name, tgt_name])
		os.system("grep -v done {0}.struct_out >> postp_out".format(outname))
		os.remove("{0}.struct_out".format(outname)) 
		os.remove("{0}.struct_out.alignment".format(outname))


if __name__ == "__main__":
    main()

