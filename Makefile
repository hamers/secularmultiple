OBJS = interface.cpp src/types.cpp src/evolve.cpp src/structure.cpp src/ODE_system.cpp src/root_finding.cpp src/newtonian.cpp src/postnewtonian.cpp src/tides.cpp src/external.cpp src/cvode/cvode.c src/cvode/cvode_dense.c src/cvode/cvode_direct.c src/cvode/cvode_io.c src/cvode/nvector_serial.c src/cvode/sundials_dense.c src/cvode/sundials_direct.c src/cvode/sundials_math.c src/cvode/sundials_nvector.c 

all: 
	$(CXX) -fPIC -shared -O1 -o libsecularmultiple.so -lm $(OBJS)
clean: 
	$(RM) libsecularmultiple.so
	
	
