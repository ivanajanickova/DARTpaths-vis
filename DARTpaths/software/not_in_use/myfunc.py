#Library of functions used in gene-pheno database management (Jens van Erp)

#given a tsv database and a list of column numbers, creates a new .tsv 
#file with the given columns in the given order.
def get_columns(fname,colnr):
	with open(fname,"r") as orig:
		destname = make_dest(fname)
		dest = open(destname,"w")
		for line in orig.readlines():
			if line[0] != "#" or line[0] != '':
				for i in range(len(colnr)):
					columns = line.split("\t")
					rightcolumn = colnr[i]
					if i+1 == len(colnr):
						dest.write(columns[rightcolumn] + "\n")
					else:
						dest.write(columns[rightcolumn] + "\t")
#specified for ensembl gene-prot databases
def get_columns2(fname):
	with open(fname,"r") as orig:
		destname = make_dest(fname)
		dest = open(destname,"w")
		for line in orig.readlines():
			if line[0] == '>':
				columns = line.split(" ")
				gene = columns[3]
				rgene = gene[5:]
				dest.write(rgene + "\t")
				prot = columns[0]
				rprot = prot[1:]
				dest.write(rprot + "\n")
#As function for be called (does not create file but returns text)
def get_columns3(text,colnr):
	linelist = []
	for line in text.split('\n'):
		if (line[0] != "#") or (line[0] != ''):
			columns = line.split("\t")
			tabs =[]
			for i in colnr:
				tabs.append(columns[i])
			newline = '\t'.join(tabs)
			linelist.append(newline)
	newtext = '\n'.join(linelist)
	return newtext

#Given a gene-database (tsv) and a ensID-to-specificID file (tsv), finds the 
#adds the ens-ID to the front of the gene-base. Very stupid algoritm 
#(nested for-loops) so can take a while
def find_genes(database,ensembl):
	with open(database,"r") as datab:
		with open(ensembl,"r") as ensem:
			destname = make_dest(database)
			dest = open(destname,"w")
			phenos = datab.readlines()
			ids = ensem.readlines()
			for an_id in ids:
				ens_id = an_id.split("\t")[0]
				own_id = an_id.split("\t")[1].rstrip()
				if own_id == "":
					continue
				for pheno in phenos:
					pheno.rstrip()
					if pheno.find(own_id) != -1:
						dest.write(ens_id + "\t" + pheno)
		
#Create a new version-filename for a destination file
def make_dest(file_name):
	if "." in file_name:
		old_name = file_name.split(".")[0]
		file_type = "." + file_name.split(".")[1]
	else:
		old_name = file_name
		file_type = ""
	version = old_name[-1]
	try:
		new_version = int(version)+1
		dest_name = old_name[:-1] + str(new_version) + file_type
	except:
		dest_name = old_name + "2" + file_type
	return dest_name

#Find all rows with (term) in a text database (db), returns rows in a list
def find_match(db,term):
	matches = []
	db.strip()
	for line in db.split('\n'):
		if line == '':
			continue
		if line[0] == "#":
			continue
		elif term in line:
			matches.append(line)
	return matches

#Give a substring (sub) within a string (full) a certain colour (col). For now, only green, yellow and red are supported.
def give_colour(full,sub,col):
	if col == "r":
		c = 1
	elif col == "g":
		c = 2
	elif col == "y":
		c = 3
	colsub = "\033[9{}m{}\033[00m".format(c,sub)
	colstr = colsub.join(full.split(sub))
	return colstr
