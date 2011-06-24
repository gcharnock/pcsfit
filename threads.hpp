

#include <boost/thread.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/function.hpp>
#include <vector>



template<typename T>
class Multithreader {
public:
    Multithreader()
		: numCPU(sysconf(_SC_NPROCESSORS_ONLN)),
		  _barrier(numCPU+1),
		  mFuncs(NULL),
          mMapTo(NULL){
		//Start the thread pool
		for(long i = 0;i<numCPU;i++) {
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
		_barrier.wait();
        //Threads do work here
		_barrier.wait();
    }
private:
	int numCPU;
    boost::barrier _barrier;

	const std::vector<boost::function<T ()> >* mFuncs;
	std::vector<T>* mMapTo;



	class WorkerThread {
	public:
		WorkerThread(long _id,Multithreader* _parent)
			: id(id),parent(_parent) {
		}
		void operator()() {
			while(true) {
				parent->_barrier.wait();
				//Do work
				for(unsigned long i = 0;i<parent->mFuncs->size();i+=parent->numCPU) {
					parent->mMapTo->at(i) = parent->mFuncs->at(i)();
				}
				parent->_barrier.wait();
			}
		}
		long id;
		Multithreader* parent;
	};

};

