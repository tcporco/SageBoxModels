# make working files in their directories.
# Note a file in a subdirectory of a working directory won't be made right,
# you'll have to construct a make -C command yourself.
wmd_files/% : /proc/uptime
	$(MAKE) -C $(dir $@) $(notdir $@)

sync :
	jekyll build --config=_config.yml,_wmd_sync.yml

# experimental pandoc pipeline
WW = /usr/local/src/workingwiki
_pandoc/Quotient.md : _pandoc/%.md : %.md.wmd
	php $(WW)/wmd/wmd.php --pre --title=Quotient --default-project-name=Box_Models --cache-dir=wmd_files --data-store=.wmd.data --modification-time=`date +%Y%m%d%H%M%s` --process-inline-math=1 --output-format=tex < $< > $@

_pandoc/Quotient.intermediate.tex : _pandoc/%.intermediate.tex : _pandoc/%.md
	pandoc -f markdown -t latex -s --listings --include-in-header=_assets/latex-header-additions.tex $< -o $@

_pandoc/Quotient.tex : _pandoc/%.tex : _pandoc/%.intermediate.tex
	php $(WW)/wmd/wmd.php --post --title=Quotient --default-project-name=Box_Models --cache-dir=wmd_files --data-store=.wmd.data --modification-time=`date +%Y%m%d%H%M%s` --output-format=tex < $< > $@

%.pdf : _pandoc/%.tex
	cd _pandoc && pdflatex $* && pdflatex $*
	mv _pandoc/$*.pdf $@

.PHONY: sync
