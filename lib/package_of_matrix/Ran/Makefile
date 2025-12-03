LIB=libran.a
OBJS=ran_imada.o
$(LIB): $(OBJS)
	ar r $(LIB) $(OBJS)
.f.o:
	$(FC) $(FLAGS) $< 
clean:
	rm $(LIB) $(OBJS)
