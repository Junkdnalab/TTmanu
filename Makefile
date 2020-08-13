clean:
	rm prereqs.R text/manuscript.html

prereqs.R :
	echo 'install.packages(c(' > prereqs.R
	grep -rih "^ *library" . | sed -E -e 's/library\(([^)]*)\).*/"\1",/' > pkglist
	cat pkglist | sed '$$s/,/))/' >> prereqs.R


text/manuscript.docx : text/manuscript.Rmd
	cd text; echo 'library(rmarkdown)\nrender("manuscript.Rmd","all")' > rcmd.R;  R CMD BATCH rcmd.R


