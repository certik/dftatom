# This part is included in all Makefiles in tests to simplify them

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) -o $@ $(OBJS) $(LDFLAGS)

test:
	./$(PROG)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<
