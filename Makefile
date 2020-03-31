# Author: Dominik Harmim <harmim6@gmail.com>

OUT := vid
TEST := test.sh
PACK := xharmi00.tgz
DOC_DIR := doc
DOC := xharmi00.pdf


.PHONY: build
build:
	@echo "Use $(TEST)." >&2
	@exit 1


.PHONY: pack
pack: $(PACK)

$(PACK):
	cp $(DOC_DIR)/$(DOC) .
	COPYFILE_DISABLE=1 tar -czf $@ $(OUT).cpp $(TEST) $(DOC)
	rm -f $(DOC)


.PHONY: clean
clean:
	rm -f $(OUT) $(PACK) $(DOC)
