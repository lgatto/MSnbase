LATEXFILES = *.aux\
        *.bbl\
        *.blg\
        *.ilg\
        *.log\
        *.nlo\
        *.nls\
        *.toc\
        *.aux\
        MSnbase-io.tex\
        MSnbase-demo.tex\
        MSnbase-development.tex\
        MSnbase-io.sh\
        MSnbase-demo.sh\
        MSnbase-development.sh\
	Rplots.pdf\
	*.dvi\
	*.map\
        *.out\
	*.figlist\
	*.dep\
	*.dpth\
	*.makefile\
	*-figure*.pdf\
	*tikzDictionary

all:
	make demo
	make devel
	make io
	make clean

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif

demo: 
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-demo.tex
	bibtex MSnbase-demo
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-demo.tex
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-demo.tex

demo2: 
	## local
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('MSnbase-demo.Rnw');"
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-demo.tex
	bibtex MSnbase-demo
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-demo.tex
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-demo.tex

io: 
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-io.tex
	bibtex MSnbase-io
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-io.tex
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-io.tex

io2: 
	## local
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('MSnbase-io.Rnw');"
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-io.tex
	bibtex MSnbase-io
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-io.tex
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-io.tex

devel: 
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-development.tex
	bibtex MSnbase-development
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-development.tex
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-development.tex

devel2: 
	## local
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('MSnbase-development.Rnw');"
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-development.tex
	bibtex MSnbase-development
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-development.tex
	"$(R_HOME)/bin/R" CMD pdflatex MSnbase-development.tex


# fake:
# 	cp ./MSnbase-demo/MSnbase-demo-fake.Rnw ./MSnbase-demo.Rnw 
# 	cp ./MSnbase-development/MSnbase-development-fake.Rnw ./MSnbase-development.Rnw 
# 	cp ./MSnbase-io/MSnbase-io-fake.Rnw ./MSnbase-io.Rnw 

# real:
# 	cp ./MSnbase-demo/MSnbase-demo.Rnw ./MSnbase-demo.Rnw 
# 	cp ./MSnbase-development/MSnbase-development.Rnw ./MSnbase-development.Rnw 
# 	cp ./MSnbase-io/MSnbase-io.Rnw ./MSnbase-io.Rnw 

clean:
	rm -f $(LATEXFILES)
	rm -f *~
	rm -f */*~
	rm -rf figure

pdfclean: 
	rm -f MSnbase-development.pdf
	rm -f MSnbase-demo.pdf
	rm -f MSnbase-io.pdf

all-clean:
	make clean
	make pdfclean
