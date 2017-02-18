TARGET   = velvob2d
MPIF90      = mpif90
CFLAGS  =  -g -fopenmp -ffree-line-length-none
OPT     = -O3

#INCLUDE = -I./include

BIN_PTH = ../bin
SRC_PTH = .
SOURCE  = $(SRC_PTH)/velvob2d.f90 $(SRC_PTH)/num2char.f $(SRC_PTH)/wavelet.f $(SRC_PTH)/p_src.f90 $(SRC_PTH)/s_src.f90 $(SRC_PTH)/getmod.f90 $(SRC_PTH)/stability.f90 $(SRC_PTH)/getfai_vsc.f90 $(SRC_PTH)/propagate_vsc_ps.f90 $(SRC_PTH)/define_pml.f90 $(SRC_PTH)/processpar.f90 $(SRC_PTH)/skiplines.f


OBJ	= $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCE)))


## Default rule executed
all: $(TARGET)
	@true
%.o: %.f90
	@echo
	@echo "=> Compiling $<"
	@$(MPIF90) $(CFLAGS) $(OPT) $(DEFINES) -c $^ -o $@ 
%.o: %.f
	@echo
	@echo "=> Compiling $<"
	@$(MPIF90) $(CFLAGS) $(OPT) $(DEFINES) -c $^ -o $@ 

$(TARGET): $(OBJ) 
	@echo
	@echo "=> Linking the target $@"
	@$(MPIF90) $(CFLAGS) $(OPT) -o $(BIN_PTH)/$@ $^ $(LIBS)
	@echo '=> Link done'
#	@echo $(OBJ) 




# the rule to clean
clean:
	@-rm *.o
	@echo "=> Clean done"

