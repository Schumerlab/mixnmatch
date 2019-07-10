#!/usr/bin/env python

# Script for extracting from a multi-fasta aligment:
# [A]ll sites (ie. do not only return variable or informative. Can still filter out N or gap-containing sites)
# [V]ariable sites  (Sites where 1 >= sites have SNPs)
# [I]nformative sites (Sites where 2 >= sites have SNPs)
# [U]ngapped alignment (remove all gap-containing sites)
# [N]-containing sites (remove sites with ambigous characters)
# [C]oordinates (return coordinates of sites rather than characters themselves)

# I overrides V

import sys, argparse
from Bio import SeqIO
from time import sleep

def Usage():
	print "Variable_sites_extraction.py -v -i -u -N -o outputfile.fasta alignment.fasta"
	sys.exit(2)

def main(argv):
	parser=argparse.ArgumentParser(description="Strip sites from multi-FASTA alignments")
	parser.add_argument("alignment")
	variable_or_informative = parser.add_mutually_exclusive_group(required=True)
	variable_or_informative.add_argument("-a", "--all", help="Return all sites. (Can still filter out N- or gap-containing sites).", action="store_true", default=False)
	variable_or_informative.add_argument("-v", "--variable", help="Return variable sites. (Sites where 1 >= isolates have SNPs)", action="store_true", default=False)
	variable_or_informative.add_argument("-i", "--informative", help="Return informative sites (Sites where 2 >= isolates have SNPs) Overrides -v", action="store_true", default=False)
	parser.add_argument("-u", "--ungapped", help="Do not return gapped sites", action="store_true", default=False)
	parser.add_argument("-N", "--N_containing", help="Do not return N-containing sites", action="store_true", default=False)
	parser.add_argument("-o", "--output_file", help="Path to output file (Default: ./output.fasta)",default="outputfile.fasta")
	parser.add_argument("-c", "--coordinates", help="Return coordinates of sites rather than characters themselves", action="store_true", default=False)
	
	args = parser.parse_args()
	
	# Needed?
	if args.variable == True:
		args.informative = False
		
	print "You have supplied the following arguments:"
	print "MSA input file:\t\t" + args.alignment
	print "Filtered MSA output file:\t\t" + args.output_file
	print "Return all sites:\t\t\t" + str(args.all)
	print "Return variable sites:\t\t\t" + str(args.variable)
	print "Return informative sites:\t\t" + str(args.informative)
	print "Filter out N-containing columns:\t" + str(args.N_containing)
	print "Filter out gapped columns:\t\t" + str(args.ungapped)
	print "Return coordinates only:\t\t" + str(args.coordinates)
					
	# FILTER ALIGNMENT:
	with open(args.alignment, "rU") as f:
		handle = SeqIO.parse(f, "fasta")
		backup = {l.id : l.seq._data for l in handle}
		dic = backup.values()
		ids = backup.keys()
		
		if not args.coordinates:
			newdic = {k : "" for k in ids}
		else:
			newdic = ""

		length = len(dic[0])
		
		print "Finished reading MSA file. Iterating over all positions and all isolates..."

		for nuc in range(length): # For each nucleotide in the alignment
			if args.all:
				include = True
			else:
				include = False
			col = "".join([x[nuc] for x in dic])
			
			# Make upper-case:
			col = col.replace("a","A").replace("c","C").replace("g","G").replace("t","T")
			
			if "N" in col and args.N_containing:
				include = False
				continue # Go to next nucleotide
			if "-" in col and args.ungapped:
				include = False
				continue # Go to next nucleotide
				
			nucleotide_counts = {}
			for char in ["A","C","G","T"]:
				charcount = col.count(char)
				if charcount > 0:
					#nucleotide_counts[char] = col.count(char)
					nucleotide_counts[char] = charcount
				
			if args.variable:
				# Then it is enough to know that more than one key exists
				if len(nucleotide_counts.keys()) > 1:
					include = True
			
			elif args.informative:
				# First we need to know that at least two different bases exists at this site
				if not len(nucleotide_counts.keys()) > 1:
					continue # Go to next nucleotide
				# Then the second most common nucleotide needs to have a value of >= 2
				if sorted(nucleotide_counts.values())[-2] >= 2:
					include = True
				if len(nucleotide_counts.keys()) == 4:
					# Sites where all bases are represented are not informative
					include = False
			
			if include:
				if not args.coordinates:
					for isolate in newdic:
						newdic[isolate] += backup[isolate][nuc]
				else:
					newdic += str(nuc) + ","
		
		# Convert to dic if coordinate mode
		if args.coordinates:
			newdic = {"C": newdic}			

		with open(args.output_file, "w") as g:
			for key in sorted(newdic):
				g.write(">" + key + "\n")
				g.write(newdic[key] + "\n")
						

if __name__ == "__main__":
	main(sys.argv[1:])
