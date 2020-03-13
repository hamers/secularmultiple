CXXSRC = interface.cpp src/types.cpp src/evolve.cpp src/structure.cpp src/ODE_system.cpp src/root_finding.cpp src/newtonian.cpp src/postnewtonian.cpp src/tides.cpp src/external.cpp src/VRR.cpp src/tools.cpp
CSRC = src/cvode/cvode.c src/cvode/cvode_dense.c src/cvode/cvode_direct.c src/cvode/cvode_io.c src/cvode/nvector_serial.c src/cvode/sundials_dense.c src/cvode/sundials_direct.c src/cvode/sundials_math.c src/cvode/sundials_nvector.c

COBJ = $(CXXSRC:.cpp=.o) $(CSRC:.c=.o)
CXXHEADERS = $(CXXSRC:.cpp=.h)
CHEADERS = $(CSRC:.c=.h)
CPPFLAGS = -fPIC -shared -O2 -Wno-comment -Wno-c++11-compat-deprecated-writable-strings

ifeq ($(DEBUG),1)
        CPPFLAGS += -DDEBUG
endif

all: $(COBJ) libsecularmultiple.so

%.o: %.c $(CHEADERS)
	@echo "Compiling C source file $< ..."
	$(CXX) -c -o $@ $< $(CPPFLAGS)

%.o: %.cpp $(CXXHEADERS)
	@echo "Compiling C++ source file $< ..."
	$(CXX) -c -o $@ $< $(CPPFLAGS)
	
libsecularmultiple.so: $(COBJ)
	@echo ""        
	@echo "Linking share library $@ ..."
	$(CXX) -o $@ $(COBJ) $(CPPFLAGS)
	@echo ""        
	@echo "The shared library $@ has been created successfully."

cleanlib: 
	$(RM) libsecularmultiple.so
clean:
	$(RM) libsecularmultiple.so *.o* src/*.o* src/cvode/*.o*
	
