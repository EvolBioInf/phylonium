include ${DISS_PATH}/Makefile

DIAGRAMS=eco29.pdf

all: manual.pdf

manual.pdf: manual.tex $(DIAGRAMS)
	pdflatex -interaction=nonstopmode $^

eco29.mat:
	phylonium ${DISS_PATH}/data/eco29/* > $@ 

eco29.tree: eco29.mat
	mat nj $^ >  $@

eco29.pdf: eco29.tree
	./etefmt.py $^ -o $@ -r eco29.names

clean: clean-latex
	setopt null_glob; $(RM) *.mat *.tree
