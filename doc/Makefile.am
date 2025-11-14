
EXTRA_DIST= dk_io.md

doc: html

html: dk_io.md
	ford -d $(abs_srcdir)/../src -o html dk_io.md

clean-local: clean-doc

clean-doc:
	rm -rf html

