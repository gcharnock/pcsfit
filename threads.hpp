
#ifndef THREADS_H
#define THREADS_H

#include <boost/thread.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/function.hpp>
#include <vector>

template<typename T>
class Multithreader {
public:
    Multithreader(bool useThreads)
		: mUseThreads(useThreads),
          numCPU(sysconf(_SC_NPROCESSORS_ONLN)),
		  _barrier(numCPU+1),
		  mFuncs(NULL) {

        mQuiting = false;

        if(!mUseThreads) {
            return;
        }

		//Start the thread pool
		for(long i = 0;i<numCPU;i++) {
            //Start up each thread and move it to the thread vector
            //(explicity copy operations are not allowed)
            mWorkers.create_thread(WorkerThread(i,this));
		}
	}

	~Multithreader() {
		//Terminate the threads
        if(!mUseThreads) {
            return;
        }
        mQuiting = true;
        _barrier.wait();
        mWorkers.join_all();
	}


    void map(const std::vector<boost::function<T ()> >& funcs,std::vector<T>& mapTo) {
		mFuncs = &funcs;
		mMapTo = &mapTo;
        mMapTo->resize(mFuncs->size());

        if(!mUseThreads) {
            //Do work
            for(unsigned long i = 0;i<mFuncs->size();i++) {
                mMapTo->at(i) = mFuncs->at(i)();
            }
        } else {
            _barrier.wait();
            //Theads to work here
            _barrier.wait();
        }
    }
private:
    bool mUseThreads;
    bool mQuiting;

	int numCPU;
    boost::barrier _barrier;

	const std::vector<boost::function<T ()> >* mFuncs;
	std::vector<T>* mMapTo;

    boost::thread_group mWorkers;

	struct WorkerThread {
		WorkerThread(long _id,Multithreader* _parent)
			: id(_id),parent(_parent) {
		}
		void operator()() {
			while(true) {
				parent->_barrier.wait();
                if(parent->mQuiting) {
                    break;
                }
				//Do work
				for(unsigned long i = id;i<parent->mFuncs->size();i+=parent->numCPU) {
					parent->mMapTo->at(i) = parent->mFuncs->at(i)();
				}
				parent->_barrier.wait();
			}
		}
		long id;
		Multithreader* parent;
	};
};

#endif
