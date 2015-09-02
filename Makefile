
CXX=c++

moltool: moltool.cpp taylor_newton.h
	$(CXX) -o $@ $<
