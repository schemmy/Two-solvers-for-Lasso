 
BOOST_ROOT=/usr/local/Cellar/boost/1.57.0/
 
project:
	g++ -O3 -std=c++11    mcx_hw1.cpp -o mcx_hw1 -L$(BOOST_ROOT)/lib -I$(BOOST_ROOT)/include -I /opt/X11/include/  -L /opt/X11/lib -lX11 -lgsl -lgslcblas -fpermissive -lboost_timer-mt -lboost_system-mt  -lboost_chrono-mt 
	./mcx_hw1	


CImgtest:
	g++ use_RGBclass.cpp -o use_RGBclass -I /opt/X11/include/ -L /opt/X11/lib -lX11 -I ~/Dropbox/PhD_Lehigh_3/407_Computational_Methods/project1/code/CImg-1.6.6_pre092915
	./use_RGBclass
