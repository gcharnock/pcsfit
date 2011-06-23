

#include <boost/thread.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/function.hpp>
#include <vector>



template<typename T>
class Multithreader {
    Multithreader()
		: numCPU(sysconf(_SC_NPROCESSORS_ONLN)),
		  _barrier(new boost::barrier(numCPU+1)),
		  mFuncs(NULL) {
		//Start the thread pool
		for(long i = 0;i<numCPU;i++) {
			boost::thread(WorkerThread(i,_barrier,*this));
		}
	}

	~Multithreader() {
		//Terminate the threads

		//todo

		delete _barrier;
	}


    void map(const std::vector<boost::function<T ()> >& funcs,std::vector<T>& mapTo) {
		mFuncs = &funcs;
		_barrier->wait();
		_barrier->wait();
    }
private:
	const std::vector<boost::function<T ()> >* mFuncs;
	std::vector<T> mMapTo;

	int numCPU;
    boost::barrier* _barrier;


	class WorkerThread {
	public:
		WorkerThread(long _id,boost::barrier* b,Multithreader _parent)
			: id(id),_barrier(b),parent(_parent) {
		}
		void operator()() {
			while(true) {
				_barrier->wait();
				//Do work
				for(long i = 0;i<parent->mFuncs.size();i+=parent->numCPU) {
					parent->mMapTo[i] = parent->mFuncs[i]();
				}
				_barrier->wait();
			}
		}
		long id;
		Multithreader& parent;
		boost::barrier* _barrier;
	};

};

