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
        MSnbase-demo.sh\
        MSnbase-development.sh\
	Rplots.pdf\
	*.dvi\
	*.map\
        *.out\

all:
	make demo
	make devel
	make clean

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif

demo: 
	mkdir -p sweave-cache/figs/
	cp ./MSnbase-demo/MSnbase-demo.Rnw . 
	##"$(R_HOME)/bin/Rscript" --vanilla -e "library(pgfSweave); pgfSweave('MSnbase-demo.Rnw');"
	"$(R_HOME)/bin/R" CMD Sweave MSnbase-demo.Rnw
	pdflatex MSnbase-demo.tex
	bibtex MSnbase-demo
	pdflatex MSnbase-demo.tex
	pdflatex MSnbase-demo.tex
	make fake

devel: 
	mkdir -p sweave-cache/figs/
	cp ./MSnbase-development/MSnbase-development.Rnw .
	##"$(R_HOME)/bin/Rscript" --vanilla -e "library(pgfSweave); pgfSweave('MSnbase-development.Rnw');"
	"$(R_HOME)/bin/R" CMD Sweave MSnbase-development.Rnw
	pdflatex MSnbase-development.tex
	bibtex MSnbase-development
	pdflatex MSnbase-development.tex
	pdflatex MSnbase-development.tex
	make fake

fake:
	cp ./MSnbase-demo/MSnbase-demo-fake.Rnw ./MSnbase-demo.Rnw 
	cp ./MSnbase-development/MSnbase-development-fake.Rnw ./MSnbase-development.Rnw 

real:
	cp ./MSnbase-demo/MSnbase-demo.Rnw ./MSnbase-demo.Rnw 
	cp ./MSnbase-development/MSnbase-development.Rnw ./MSnbase-development.Rnw 

clean:
	rm -f $(LATEXFILES)
	rm -rf sweave-cache
	rm -f *~

pdfclean: 
	rm -f MSnbase-development.pdf
	rm -f MSnbase-demo.pdf

all-clean:
	make clean
	make pdfclean
