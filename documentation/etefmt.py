#!/usr/bin/env python

from ete3 import Tree, TreeStyle, NodeStyle, TextFace, AttrFace, faces
from pathlib import Path
import argparse
import csv
import os
import sys

out_stem = "etefmt" # not actual default
out_format = "pdf" # not actual default
names = {}

def layout(node):
	if node.is_leaf(): 
		pretty_name = names[node.name] if node.name in names else node.name
		# derp = AttrFace("name", fsize=10, ftype="Arial", fgcolor="black")
		pretty_name_face = faces.TextFace(pretty_name)
		pretty_name_face.margin_left = 5
		faces.add_face_to_node(pretty_name_face, node, column=0)

def read_tree(file_name):
	if file_name == "-":
		file_name = "/dev/stdin"

	t = Tree(file_name);
	
	# ignore any errors with midpoint rooting
	R = t.get_midpoint_outgroup()
	try:
		t.set_outgroup(R)
	except Exception as e:
		pass

	t.ladderize(direction=1)
	return t
	

def main(file_name):
	t = read_tree(file_name)

	nstyle = NodeStyle()
	nstyle["size"] = 0

	ts = TreeStyle()
	ts.show_leaf_name = False
	ts.scale = 25000 # pixels per branch length unit
	ts.layout_fn = layout

	for n in t.traverse():
		n.set_style(nstyle)

	# t.show(tree_style=ts)
	t.render(out_stem + "." + out_format, tree_style=ts)

def read_names(names_file):
	global names
	with open(names_file) as csvfile:
		namesreader = csv.reader(csvfile, delimiter=',')
		for row in namesreader:
			names[row[0]] = row[1]
	return names

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Render a phylogeny to an image.')
	parser.add_argument('file', nargs='?', help='Name of file with phylogeny')
	parser.add_argument('-o', dest='image_file', help='Image file name')
	parser.add_argument('-r', dest='names_file', help='File with ID and real name (CSV)')
	args = parser.parse_args()

	if args.names_file:
		read_names(args.names_file)

	# parse output file name
	if args.image_file:
		# split arg
		out_stem = Path(args.image_file).stem
		out_format = Path(args.image_file).suffix[1:]
	else:
		# split input file name
		print(args.file)
		inputfile = Path(args.file if args.file else "etefmt.pdf")
		out_stem = Path(inputfile).stem
		out_format = "pdf"

	if args.file:
		main(args.file)
	elif not os.isatty(0):
		main("-")
	else:
		parser.print_help()
		


