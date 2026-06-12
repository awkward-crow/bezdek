PKG     := bezdek
VERSION := 0.1.0
TARBALL := $(PKG)_$(VERSION).tar.gz

.PHONY: build check install clean

build:
	R CMD build .

check: build
	R CMD check $(TARBALL)

install:
	R CMD INSTALL .

clean:
	$(RM) $(TARBALL)
	$(RM) -r $(PKG).Rcheck
	$(RM) src/*.o src/*.so src/*.dll
