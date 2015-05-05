# make working files in their directories.
# Note a file in a subdirectory of a working directory won't be made right,
# you'll have to construct a make -C command yourself.
wmd_files/% : /proc/uptime
	$(MAKE) -C $(dir $@) $(notdir $@)

sync :
	jekyll build --config=_config.yml,_wmd_sync.yml

# experimental pandoc pipeline
# TODO: need to get the project name and title right from the YAML
# ideally the rest of the YAML too, though we can parasitize it from
# .wmd.data if we don't mind being kind of low-quality (this requires
# .wmd.data to be left over from a jekyll build operation, and all the
# source files from other pages to be freshly synced, because jekyll
# and pandoc disagree on locations within the page text, with and without
# the YAML header included).
WW = /usr/local/src/workingwiki
PROJECT=Notes
TITLE="Definitions of Box Model objects"
_pandoc/%.md : %.md.wmd wmd_files/.workingwiki/.wmd.data
	php $(WW)/wmd/wmd.php --pre --title=$(TITLE) --default-project-name=$(PROJECT) --cache-dir=wmd_files --data-store=.wmd.data --modification-time=`date +%Y%m%d%H%M%s` --process-inline-math=1 --output-format=tex < $< > $@

_pandoc/%.intermediate.tex : _pandoc/%.md
	pandoc -f markdown -t latex -s --listings --include-in-header=_assets/latex-header-additions.tex $< -o $@

_pandoc/%.tex : _pandoc/%.intermediate.tex
	php $(WW)/wmd/wmd.php --post --title=$(TITLE) --default-project-name=$(PROJECT) --cache-dir=wmd_files --data-store=.wmd.data --persistent-data-store --modification-time=`date +%Y%m%d%H%M%s` --output-format=tex < $< > $@

%.pdf : _pandoc/%.tex
	cd _pandoc && pdflatex $* && pdflatex $*
	mv _pandoc/$*.pdf $@

wmd_files/.workingwiki/.wmd.data : *.md.wmd # */*.md.wmd
	$(MAKE) sync

.PRECIOUS: _pandoc/%.md _pandoc/%.tex

.PHONY: sync
