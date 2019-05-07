all:
	$(MAKE) build
	go install
build:
	mkdir _build/ libs/
	cmake -B _build/ -S ..
	$(MAKE) -C _build/
	cp _build/libsymspg.a libs/
	rm -rf _build

clean:
	rm -rf _build/
	rm -rf libs/
