MAIN = mailler

# Le programme principal doit s'appeler $(MAIN).f90, les modules
# doivent s'appeler mod_XXX.f90 et l'executable final
# s'appelera $(MAIN) en majuscules.
SORTIE = $(shell echo $(MAIN) | tr a-z A-Z)
COMP = gfortran
OPTS = -O5 -fbackslash -fintrinsic-modules-path /usr/local/lib/mods/
OBJS_f90 = $(filter-out calcul_taux.f90 mod_interface_hsl.f90,$(wildcard *.f90))
OBJS = $(patsubst %.f90,%.o,$(OBJS_f90))
LIBS = 

all: $(SORTIE)

##
# BINAIRES RELOGEABLES
# Pour les modules, les .mod ne sont jamais generes eux-memes, mais
# comme ils dependent des .o, la creation de ces derniers mettra a jour
# les .mod correspondants
#mod_%.mod: mod_%.o ;
#mod_%.o: mod_%.f90
#	$(COMP) $(OPTS) -c $<
#mod_%.o mod_%.mod: mod_%.f90
#	$(COMP) $(OPTS) -c $<
# We 'touch' the .mod because gfortran does not re-create
# it if it hasn't changed, which contradicts 'make' view
# of the world (cmake doesn't seem to have this issue)
mod_maillage.o mod_maillage.mod: mod_maillage.f90 mod_variables.mod mod_dessin.mod;
	$(COMP) $(OPTS) -c $<
	touch mod_maillage.mod
mod_dessin.o mod_dessin.mod: mod_dessin.f90 mod_variables.mod;
	$(COMP) $(OPTS) -c $<
	touch mod_dessin.mod
mod_variables.o mod_variables.mod: mod_variables.f90 ;
	$(COMP) $(OPTS) -c $<
	touch mod_variables.mod
mod_castests.o mod_castests.mod: mod_castests.f90 mod_variables.mod;
	$(COMP) $(OPTS) -c $<
	touch mod_castests.mod
mod_reslin.o mod_reslin.mod: mod_reslin.f90
	$(COMP) $(OPTS) -c $<
	touch mod_reslin.mod
mod_interface_hsl.o mod_interface_hsl.mod: mod_interface_hsl.f90
	$(COMP) $(OPTS) -c $<
	touch mod_interface_hsl.mod


# Le .o du programme principal depend de son .f90 et des mod_XXX.mod
# (seuls necessaires a creer le binaire relogeable principal)
$(MAIN).o: $(MAIN).f90 $(patsubst mod_%.f90,mod_%.mod,$(wildcard mod_*.f90))
	$(COMP) $(OPTS) -c $<

##
# EDITION DES LIENS
# La generation de l'executable depend de tous les .o (modules et programme principal)
$(SORTIE): $(OBJS)
	$(COMP) $(OPTS) $(OBJS) $(LIBS) -o $@

##
# EXECUTION
.PHONY: run
run: $(SORTIE)
	./$(SORTIE)

##
# Calcul des taux de convergence
calcul_taux: CALCUL_TAUX
CALCUL_TAUX: calcul_taux.f90
	$(COMP) $(OPTS) $< -o CALCUL_TAUX

calcul_taux_run: calcul_taux
	./CALCUL_TAUX

##
# NETTOYAGE
.PHONY : clean
clean: 
	@echo "Nettoyage";
	rm -f *.o *.mod *.fig *.vtk;
  
