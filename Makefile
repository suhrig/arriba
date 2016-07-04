.PHONY: all
all: ariba extract_read-through_fusions

.PHONY: ariba
ariba:
	$(MAKE) -C source ariba
	cp source/ariba .

.PHONY: extract_read-through_fusions
extract_read-through_fusions:
	$(MAKE) -C source extract_read-through_fusions
	cp source/extract_read-through_fusions .

.PHONY: clean
clean:
	$(MAKE) -C source clean
	$(MAKE) -C samtools-1.3 clean
	$(MAKE) -C samtools-1.3/htslib-1.3 clean

