
# to make a single target, you can use make wmd_files/Project/target
# or more usefully, make sync wmd_files/Project/target
# using the rule below.
# that puts the output of make on stdout.  To get the output in the
# .make.log file, do make wmd_files/Project/target.make.log
# using this rule.  Note sync is redundant because this includes the
# sync operation.
%.make.log : /proc/uptime
	php $(WW_DIR)/wmd/wmd.php --post --cache-dir=wmd_files --default-project-name=$(subst /,,$(subst wmd_files/,,$(dir $*))) --make-single-file=$(notdir $*)

kdjlfj:
	echo wmd_make_single_file: $(notdir $*) > _make_single_file.yml
	jekyll build --config=_config.yml,_make_single_file.yml
	$(RM) _make_single_file.yml

# make working files in their directories.
# Note a file in a subdirectory of a working directory won't be made right,
# you'll have to construct a make -C command yourself.
wmd_files/% : /proc/uptime
	$(MAKE) -C $(dir $@) $(notdir $@)

sync :
	jekyll build --config=_config.yml,_wmd_sync.yml

# generate CSS from WW CSS
WW_DIR = /usr/local/src/workingwiki

WW_CSS_TO_USE = $(patsubst %,$(WW_DIR)/resources/%,ext.workingwiki.latexml.css ext.workingwiki.latexml.customization.css)

css/auto-generated-from-ww.css : $(WW_CSS_TO_USE)
	cat $(WW_CSS_TO_USE) > $@

# to rebuild a single page, i.e. rebuild the site while only doing make
# operations on a single page
# $ make _site/page.html
_site/%.html : /proc/uptime
	echo wmd_make_page: $*.md > _make_page.yml
	jekyll build --config=_config.yml,_make_page.yml
	$(RM) _make_page.yml

# experimental pandoc pipeline

WW = /usr/local/src/workingwiki
_pandoc/%.md : %.md.wmd wmd_files/.workingwiki/.wmd.data
	php $(WW)/wmd/wmd.php --pre --title='$(TITLE)' --default-project-name=$(PROJECT) --cache-dir=wmd_files --data-store=.wmd.data --modification-time=`date +%Y%m%d%H%M%s` --process-inline-math=1 --output-format=tex < $< > $@

_pandoc/%.intermediate.tex : _pandoc/%.md
	pandoc -f markdown -t latex -s -S --listings --include-in-header=_assets/latex-header-additions.tex --filter pandoc-citeproc $< -o $@

_pandoc/%.tex : _pandoc/%.intermediate.tex
	php $(WW)/wmd/wmd.php --post --title='$(TITLE)' --default-project-name=$(PROJECT) --cache-dir=wmd_files --data-store=.wmd.data --persistent-data-store --modification-time=`date +%Y%m%d%H%M%s` --output-format=tex < $< > $@

_pandoc/Definitions.intermediate.tex : _pandoc/%.intermediate.tex : _pandoc/%.md box.bib
	cp box.bib _pandoc
	pandoc -f markdown -t latex -s -S --listings --include-in-header=_assets/latex-header-additions.tex --filter pandoc-citeproc $< -o $@
	$(RM) _pandoc/box.bib

#Definitions.pdf : %.pdf : _pandoc/%.tex box.bib
#	cd _pandoc && pdflatex $* && bibtex $* && pdflatex $* && pdflatex $*
#	mv _pandoc/$*.pdf $@

%.pdf : _pandoc/%.tex
	cd _pandoc && pdflatex $* && pdflatex $*
	mv _pandoc/$*.pdf $@

wmd_files/.workingwiki/.wmd.data : *.md.wmd # */*.md.wmd
	$(MAKE) sync

# TODO: need to get the project name and title right from the YAML
# ideally the rest of the YAML too, though we can parasitize it from
# .wmd.data if we don't mind being kind of low-quality (this requires
# .wmd.data to be left over from a jekyll build operation, and all the
# source files from other pages to be freshly synced, because jekyll
# and pandoc disagree on locations within the page text, with and without
# the YAML header included).
_pandoc/Measles.% : PROJECT=Measles 
_pandoc/Measles.% : TITLE="Subcritical Measles Outbreak Size"
_pandoc/Definitions.% : PROJECT=Notes
_pandoc/Definitions.% : TITLE="Definitions of Box Model objects"
_pandoc/sde.% : PROJECT=SDE
_pandoc/sde.% : TITLE=SDE
_pandoc/LeprosyDiagram.% : PROJECT=Leprosy
_pandoc/LeprosyDiagram.% : TITLE="Leprosy Diagram"

.PRECIOUS: _pandoc/%.md _pandoc/%.tex

.PHONY: sync
