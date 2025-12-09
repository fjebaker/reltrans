BUILD  := build
ROOTDIR := ..
HEADAS_LIB := ${HEADAS}/lib
HEADAS_INCLUDE := ${HEADAS}/include

FC := gfortran

CFLAGS := -g -fno-omit-frame-pointer

FFLAGS := -DHAVE_INLINE -g -fPIC -fno-automatic -fno-second-underscore \
		  -fopenmp \
		  -O3 \
		  -I$(BUILD)/include \
		  -I$(HEADAS_INCLUDE) \
		  -J$(BUILD)/cache \
		  -I$(BUILD)/cache

LDFLAGS := -lXSFunctions -lXSModel -lfftw3 -lm -lcfitsio -lpthread \
		   -L$(BUILD)/lib -L$(HEADAS_LIB) -shared -export-dynamic

all: $(BUILD) $(BUILD)/lib/libreltrans.so

exe: $(BUILD)/bin/relcli

$(BUILD)/bin/relcli: ./cli.c $(BUILD)/lib/libreltrans.so
	$(CC) $(CFLAGS) cli.c -o $(BUILD)/bin/relcli -Wl,-rpath,$(BUILD)/lib \
		-L$(BUILD)/lib -lgfortran -lc -lreltrans -lm -lmvec

$(BUILD)/lib/libreltrans.so: $(BUILD)/cache/wrappers.o
	$(LD) $(LDFLAGS) $(BUILD)/cache/wrappers.o -o $@

$(BUILD)/cache/%.o: $(ROOTDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)/bin $(BUILD)/lib $(BUILD)/include $(BUILD)/share $(BUILD)/cache

.PHONY: clean
clean:
	rm -rf $(BUILD)

.PHONY: format
format:
	clang-format -i ./cli.c
