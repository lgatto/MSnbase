all:
	make MSnbase-io-in.png
	make MSnbase-io-out.png

MSnbase-io-in.pdf: MSnbase-io-in.tex
	pdflatex MSnbase-io-in.tex

MSnbase-io-in.png: MSnbase-io-in.pdf
	convert -density 150 MSnbase-io-in.pdf MSnbase-io-in.png

MSnbase-io-out.pdf: MSnbase-io-out.tex
	pdflatex MSnbase-io-out.tex

MSnbase-io-out.png: MSnbase-io-out.pdf
	convert -density 150 MSnbase-io-out.pdf MSnbase-io-out.png

clean:
	rm -f *~
	rm -f MSnbase-io-in.pdf MSnbase-io-out.pdf
	rm -f *aux *log
