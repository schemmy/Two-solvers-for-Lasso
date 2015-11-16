 
BOOST_ROOT=/usr/local/Cellar/boost/1.57.0
 
project:
	g++ -O3 -std=c++11    main.cpp -o main -L$(BOOST_ROOT)/lib -I$(BOOST_ROOT)/include \
		-I /opt/X11/include/  -L /opt/X11/lib -lX11 -lgsl -lgslcblas -fpermissive \
		 -lboost_timer-mt -lboost_system-mt  -lboost_chrono-mt 
	./main -A photos/test4 -M 2
