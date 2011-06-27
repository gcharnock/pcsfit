
#include <boost/version.hpp>
#include <boost/thread.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/function.hpp>
#include <iostream>
#include <vector>



template<typename T>
class Multithreader {
public:
    Multithreader()
		: numCPU(sysconf(_SC_NPROCESSORS_ONLN)),
		  _barrier(numCPU+1),
		  mFuncs(NULL),
		  mMapTo(NULL) {
		std::cout << "BOOST_VERSION = " << BOOST_VERSION << std::endl;
		//Start the thread pool
		for(long i = 0;i<numCPU;i++) {
			std::cout << "Starting thread "<<i<< std::endl;
			boost::thread(WorkerThread(i,this));
		}
	}

	~Multithreader() {
		//Terminate the threads

		//todo
	}


    void map(const std::vector<boost::function<T ()> >& funcs,std::vector<T>& mapTo) {
		mFuncs = &funcs;
		mMapTo = &mapTo;		

		std::cout << "Threads have work..." << std::endl;
		_barrier.wait();
		//Workers are working at this point
		std::cout << "Threads are working..." << std::endl;
		_barrier.wait();
		std::cout << "Threads have finished" << std::endl;
    }
private:
	const std::vector<boost::function<T ()> >* mFuncs;
	std::vector<T>* mMapTo;
    boost::barrier _barrier;

	size_t numCPU;

	class WorkerThread {
	public:
		WorkerThread(long _id,Multithreader* _parent)
			: id(_id),parent(_parent) {
		}
		void operator()() {
			while(true) {
				std::cout << "Thread "<<id<<" waiting for work" << std::endl;
				parent->_barrier.wait();
				std::cout << "Thread "<<id<<" has work" << std::endl;
				//Do work
				for(long i = id;i<parent->mFuncs->size();i+=parent->numCPU) {
					parent->mMapTo->at(i) = parent->mFuncs->at(i)();
				}
				parent->_barrier.wait();
				std::cout << "Thread "<<id<<" finished the work" << std::endl;
			}
		}
		long id;
		Multithreader* parent;
	};

};

