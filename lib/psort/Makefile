.PHONY: clean

PROGRAM = libpsort.a

LDFLAGS = -I../ohhelp/ohhelp-1.1.1/
CFLAGS := $(LDFLAGS) $(CFLAGS)

OBJS = psort-struct.o


$(PROGRAM): $(OBJS)
	ar rc $@ $^

$(OBJS):%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:;
	rm $(OBJS)
	rm $(PROGRAM)