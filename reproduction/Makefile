R_OPTS=--no-save --no-restore --no-init-file --no-site-file

reproduction.html: reproduction.Rmd
	R $(R_OPTS) -e "rmarkdown::render('$<')"
