#######################################################
#
# Welcome to the Makefile!
#  > if you are looking for compile time options
#  > (including compiler definition)
#  > you do not need to modify anything in this file
#  > check 'options.make' file instead
#
#######################################################

 PROGRAMNAME=patchwork

 FC=ifort
 SRCDIR=src
 MODDIR=mod
 OBJDIR=obj
# MODF=$(SRCDIR)/vars.ftn

# SRC=$(sort $(wildcard ls -1 $(SRCDIR)/*f90))
 SRC :=     math.f90 \
            stretching.f90 \
            io.f90 \
            data_types.f90 \
            cubic_mod.f90 \
	        main.f90 \

 OBJ=$(patsubst %.f90,$(OBJDIR)/%.o,$(notdir $(SRC)))

# Intel Fortran Compiler flags
FLAGS= -traceback #-g -O0 -warn all -debug extended -nofpp -fpe0 -double-size 64 -traceback -debug -check 
SETMOD=-module mod
NOFX=-nofixed

# gfortran flags
# FLAGS= -ffree-line-length-200

###################################################################
######################      RECIPES       #########################
###################################################################


$(PROGRAMNAME): $(OBJDIR) $(MODDIR) $(OBJ)
	@echo " "
	@echo "linking..."
	@echo " "
	$(FC) $(FLAGS) $(SETMOD) -o $@ $(OBJ) 
	@echo " "
	@echo "all done!"
	@echo " "

$(OBJDIR):
	mkdir $(OBJDIR)

$(MODDIR):
	mkdir $(MODDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.f90  
	$(FC) $(FLAGS) $(SETMOD) -c $< -o $@ 2>errors

#------------------------------------------------------------------
.PHONY: clean tools

clean:
	@rm -vf $(PROGRAMNAME)
	@rm -vf $(MODDIR)/*genmod*
	@rm -vf $(MODDIR)/*.mod
	@rm -vf $(OBJDIR)/*.o
	rmdir  $(MODDIR)
	rmdir  $(OBJDIR)

