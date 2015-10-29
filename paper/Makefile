SOURCES=$(wildcard *.tex)

all: paper.pdf

paper.pdf: benchmark.pdf comparison.pdf merge_example.pdf paper.bbl $(SOURCES)
	pdflatex paper
	pdflatex paper

bib:
	pdflatex paper
	bibtex paper

paper.bbl: paper.bib
	pdflatex paper
	bibtex paper

benchmark.pdf: benchmark.csv benchmark.R
	R --slave --args benchmark < benchmark.R

comparison.pdf: comparison.csv comparison.R
	R --slave --args comparison < comparison.R

merge_example.pdf: merge_example.ipe
	ipetoipe -pdf merge_example.ipe

clean:
	rm -f *.aux *.bak *.bbl *.blg *.log *.pdf
