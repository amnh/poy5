SVNVERSION := svnversion

QuickStart.pdf: QuickStart.tex
	pdflatex $<; pdflatex $<; pdflatex $< 

commands.pdf: commands.tex allcommands.tex poy3commands.tex  QuickStart.tex poylibrary.bib poyheuristics.tex
	echo '\newcommand{\buildnumber}{' `$(SVNVERSION) -n ./` " }" | sed "s/ //g" > version.tex; pdflatex commands.tex ; bibtex commands; pdflatex commands.tex ; pdflatex commands.tex ; makeindex general ; makeindex poy3 ; pdflatex commands.tex

clean:
	rm -f *.aux *.pdf *.log *.idx *.ind *.out *.toc *.ilg *.bbl *.blg

