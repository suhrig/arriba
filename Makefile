.PHONY: all
all: arriba extract_reads

.PHONY: arriba
arriba:
	$(MAKE) -C source arriba
	cp source/arriba .

.PHONY: extract_reads
extract_reads:
	$(MAKE) -C source extract_reads
	cp source/extract_reads .

.PHONY: clean
clean:
	$(MAKE) -C source clean
	$(MAKE) -C samtools-1.3 clean
	$(MAKE) -C samtools-1.3/htslib-1.3 clean

