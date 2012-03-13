#!/usr/bin/env python
import re, os


file = open("bc-30.db")
pattern = r'^(name:)\s(\w+)'
db_tail = 'db'

for line in file:
	match = re.search(pattern, line)
	if match:
		try:
			file_out
		except NameError:
			print "first line"
		else:
			file_out.close()
		print match.group(2)[0:4]
		file_out_name = "dbis/" + match.group(2) + '.db'
		file_out = open(file_out_name, 'w')
		command = "./pdbdownload.pl " + match.group(2)[0:4]
		os.system(command)
	file_out.write(line)
file.close()
file_out.close()


