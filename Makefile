# make working files in their directories.
# Note a file in a subdirectory of a working directory won't be made right,
# you'll have to construct a make -C command yourself.
wmd_files/% : /proc/uptime
	$(MAKE) -C $(dir $@) $(notdir $@)

sync :
	jekyll build --config=_config.yml,_wmd_sync.yml
