.PHONY: all
all: arriba extract_read-through_fusions

.PHONY: arriba
arriba:
	$(MAKE) -C source arriba
	cp source/arriba .

.PHONY: extract_read-through_fusions
extract_read-through_fusions:
	$(MAKE) -C source extract_read-through_fusions
	cp source/extract_read-through_fusions .

.PHONY: clean
clean:
	$(MAKE) -C source clean
	$(MAKE) -C samtools-1.3 clean
	$(MAKE) -C samtools-1.3/htslib-1.3 clean

