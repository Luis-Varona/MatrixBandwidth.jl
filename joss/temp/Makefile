all: paper.pdf

paper.pdf: ../paper.md ../paper.bib latex.template
	cd ../ && pandoc --citeproc --bibliography paper.bib paper.md \
	--template temp/latex.template -o temp/paper.pdf \
	--pdf-engine=lualatex

clean:
	rm -f paper.pdf

.PHONY: clean
