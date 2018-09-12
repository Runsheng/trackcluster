CC = gcc
CFLAGS := -Wall -O3 -pipe #-pg
LIB = libssw.so

.PHONY: all clean


$(LIB): ssw.c ssw.h
	$(CC) $(CFLAGS) -fPIC -shared -rdynamic -o $@ $<

clean:
	-rm -f $(LIB) 
