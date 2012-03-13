#!/usr/bin/env python
import os
import sys

def inOrder(line, order):
	in_order = True
	words = line.split()
	key = words[5]
	value = int(words[1][0:-1].split("_")[1])
	
	if key in order:
		if value < order[key]:
			in_order = False
		else :
			order[key] = value
	else :
		order[key] = value
	return in_order

def chaseOrder(file):
	order = {}
	for line in file:
		if 'select qry' in line:
			if inOrder(line, order) == False:
				return False
			
		elif 'select tgt' in line:
			if inOrder(line, order) == False:
				return False
	return True

def main():                                        

# parsing input values
                         
	basedir = sys.argv[1]
	for root, dirs, files in os.walk(basedir):
		os.chdir(root)
		print len(files)
		for file in files:
			if '.pml' in file:
				f = open(file, "r")
				if chaseOrder(f) == False:
					print file
				
				


if __name__ == "__main__":
    main()

	