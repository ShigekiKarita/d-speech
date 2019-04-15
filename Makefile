.PHONY: doc

adrdox:
	git clone https://github.com/adamdruppe/adrdox.git --depth 1

adrdox/doc2: adrdox
	cd adrdox; make

doc: adrdox/doc2
	./adrdox/doc2 -u -i source

