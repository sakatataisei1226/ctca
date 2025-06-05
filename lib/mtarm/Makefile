.PHONY: clean

PROGRAM = libmtarm.a

OBJS = mt19937arm.o

CFLAGS := $(CFLAGS)

$(PROGRAM): $(OBJS)
	ar rc $@ $^

$(OBJS):%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:;
	rm $(OBJS)
	rm $(PROGRAM)
