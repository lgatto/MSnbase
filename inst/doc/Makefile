LATEXFILES = *.aux\
        *.bbl\
        *.blg\
        *.ilg\
        *.log\
        *.nlo\
        *.nls\
        *.toc\
        *.aux\
        MSnbase-demo.tex\
        MSnbase-development.tex\
	*.map\
        *.out

make all:
	make all-clean
	make demo
	make devel

demo: 
	R CMD pgfsweave --graphics-only MSnbase-demo.Rnw
	pdflatex MSnbase-demo.tex
	bibtex MSnbase-demo
	pdflatex MSnbase-demo.tex
	pdflatex MSnbase-demo.tex

devel: 
	R CMD pgfsweave --graphics-only MSnbase-development.Rnw
	pdflatex MSnbase-development.tex
	bibtex MSnbase-development
	pdflatex MSnbase-development.tex
	pdflatex MSnbase-development.tex

clean:
	rm -f $(LATEXFILES)
	rm -f *~
	rm -f MSnbase-demo.sh
	rm -f MSnbase-development.sh
	rm -f Rplots.pdf

all-clean:
	make clean
	rm -f MSnbase-development.pdf
	rm -f MSnbase-development.dvi
	rm -f MSnbase-demo.pdf
	rm -f MSnbase-demo.dvi
	rm -f sweave-cache/figs/*

