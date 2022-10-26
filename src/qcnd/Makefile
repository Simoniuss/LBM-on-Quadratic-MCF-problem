CFLAGS = -O -DDIMACS

TARGETS = pargen  qfcgen

.PHONY: all clean

%: %.c
	$(CXX) -o $@ $<

all: $(TARGETS)
	$(MAKE) -C ./netgen

clean : 
	rm -f $(TARGETS)
	$(MAKE) -C ./netgen clean

