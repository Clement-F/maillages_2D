MAIN = main

# Le programme principal doit s'appeler $(MAIN).F90, les modules
# doivent s'appeler mod_XXX.F90 et l'executable final
# s'appelera $(MAIN) en majuscules.
SORTIE = $(shell echo $(MAIN) | tr a-z A-Z)
COMP = gfortran
OPTS = -O5 -fbackslash -fintrinsic-modules-path /usr/local/lib/mods/
OBJS_f90 = $(filter-out calcul_taux.F90 mod_interface_hsl.F90,$(wildcard *.F90))
OBJS = $(patsubst %.F90,%.o,$(OBJS_f90))
LIBS = 


FFLAGS  =  -fcheck=all -g -C -fbacktrace -Wtarget-lifetime -Wcompare-reals -Wrealloc-lhs-all \
			-Wrealloc-lhs -Wintrinsic-shadow -Wunderflow -Wintrinsics-std \
			-Wimplicit-procedure -Wconversion -Wconversion-extra -Wline-truncation \
			-Wcharacter-truncation -Wc-binding-type -Waliasing -ffpe-trap=zero,overflow,underflow
LDFLAGS = -g -C -fbacktrace -Wtarget-lifetime -Wcompare-reals -Wrealloc-lhs-all \
			-Wrealloc-lhs -Wintrinsic-shadow -Wunderflow -Wintrinsics-std \
			-Wimplicit-procedure -Wconversion -Wconversion-extra -Wline-truncation \
			-Wcharacter-truncation -Wc-binding-type -Waliasing -ffpe-trap=zero,overflow,underflow 

all: $(SORTIE)

##
# BINAIRES RELOGEABLES
# Pour les modules, les .mod ne sont jamais generes eux-memes, mais
# comme ils dependent des .o, la creation de ces derniers mettra a jour
# les .mod correspondants
#mod_%.mod: mod_%.o ;
#mod_%.o: mod_%.F90
#	$(COMP) $(FFLAGS) $(OTPS) -c $<
#mod_%.o mod_%.mod: mod_%.F90
#	$(COMP) $(FFLAGS) $(OTPS) -c $<
# We 'touch' the .mod because gfortran does not re-create
# it if it hasn't changed, which contradicts 'make' view
# of the world (cmake doesn't seem to have this issue)

mod_mesh.o mod_mesh.mod: mod_mesh.F90 mod_var.mod mod_divers.mod;
	$(COMP) $(FFLAGS) $(OTPS) -c $<
	touch mod_mesh.mod

mod_divers.o mod_divers.mod: mod_divers.F90 precis.mod mod_var.mod;
	$(COMP) $(FFLAGS) $(OTPS) -c $<
	touch mod_divers.mod

precis.o precis.mod: precis.F90;
	$(COMP) $(FFLAGS) $(OTPS) -c $<

mod_var.o mod_var.mod: mod_var.F90 precis.mod ;
	$(COMP) $(FFLAGS) $(OTPS) -c $<
	touch mod_var.mod

# Le .o du programme principal depend de son .F90 et des mod_XXX.mod
# (seuls necessaires a creer le binaire relogeable principal)
$(MAIN).o: $(MAIN).F90 $(patsubst mod_%.F90,mod_%.mod,$(wildcard mod_*.F90))
	$(COMP) $(FFLAGS) $(OTPS) -c $<

##
# EDITION DES LIENS
# La generation de l'executable depend de tous les .o (modules et programme principal)
$(SORTIE): $(OBJS)
	$(COMP) $(FFLAGS) $(OTPS) $(OBJS) $(LIBS) -o $@

##
# EXECUTION
.PHONY: run
run: $(SORTIE)
	./$(SORTIE)

##
# NETTOYAGE
.PHONY : clean
clean: 
	@echo "Nettoyage";
	rm -f *.o *.mod *.fig *.vtk;
  
