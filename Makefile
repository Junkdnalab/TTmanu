clean:
	rm prereqs.R text/manuscript.html

prereqs.R :
	echo 'install.packages(c(' > prereqs.R
	grep -rih "^library" . | sed -E -e 's/library\(([^)]*)\).*/"\1",/' > pkglist
	cat pkglist | sed '$$s/,/))/' >> prereqs.R


text/manuscript.html : text/manuscript.Rmd
	echo 'library(rmarkdown)\nrender("text/manuscript.Rmd","html_document")' > rcmd.R;  R CMD BATCH rcmd.R


