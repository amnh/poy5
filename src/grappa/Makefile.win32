all: $(DEPEND) $(FILES)

COBJS = uf.o \
	invdist.o correction.o binencode.o lists.o med_util.o \
	all_sorting_reversals.o vertex_factory.o hashtable.o \
	priority_stack.o sorting_reversal_median.o \
	inversion_median.o inversion_median_alberto.o grappa_interface.o 

libcamlgrappa.a: $(COBJS)
	rm -f $@
	ar rc $@ $(COBJS)

test_grappa: grappa.cmo	test_grappa.ml
	$(OCAMLC) $(CCOPT) $(OCAML_INCLUDES) -custom -cclib -g -cclib -L./ \
	-cclib -lcamlgrappa grappa.cma test_grappa.ml -o $@

.SUFFIXES: .ml .mli .cmo .cmx .cmi

$(DEPEND):	*.ml *.mli
	ocamldep -I $(INCLUDEOCAMLMPI) *.ml *.mli > $(DEPEND)

.ml.cmo:
	$(OCAMLC) -g -c $<

%.cmx: %.ml
	$(OCAMLOPT) -c $(CCOPT) $< $(PROFILE)

.ml.cmx:
	$(OCAMLC) -c $(CCOPT) $< $(PROFILE)

.mli.cmi:
	$(OCAMLC) -c -I $(INCLUDE) $(DEBUG) $<

.c.o:
	$(OCAMLOPT) $(CCOPT) -c $< $(PROFILE)


grappa.cmo: grappa.cmi

grappa.cmx: grappa.cmi

grappa.cma:  libcamlgrappa.a grappa.cmo 
	$(OCAMLC) -a -o grappa.cma -custom -cclib -L./ -cclib -lcamlgrappa grappa.cmo

grappa.cmxa:  libcamlgrappa.a grappa.cmx
	$(OCAMLOPT) -a -o grappa.cmxa -cclib -L./ -cclib -lcamlgrappa grappa.cmx
	
clean:	
	rm -f *.cmi *.cmo *.cmx *.o test_grappa *.cma *.cmxa *.a *.annot
