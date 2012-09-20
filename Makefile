SRC_DIR = src

SRC_C_FILES = src/constraints.c src/RNAmutants.c src/util.c src/sampling.c src/mfe_backtrack.c src/energy_functions.c src/reader_energy_params.c src/hybrid.c

SRC_H_FILES = src/util.h src/constraints.h src/sampling.h src/reader_energy_params.h src/energy_params_tables.h src/mfe_backtrack.h src/RNAmutants.h src/energy_functions.h src/hybrid.h


all	: RNAmutants

RNAmutants	: $(SRC_C_FILES) $(SRC_H_FILES) 
	(cd $(SRC_DIR)  &&  $(MAKE) $@) ||  exit 1;\
	(mv $(SRC_DIR)/RNAmutants  .) ||  exit 1

clean	: 
	(cd $(SRC_DIR)  &&  $(MAKE) $@)  ||  exit 1;\
	(rm -f *~) || exit 1
