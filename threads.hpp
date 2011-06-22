

#include <boost/thread.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/function.hpp>
#include <vector>

template<typename T>
class Multithreader {
    Multithreader() {
        int numCPU = sysconf(_SC_NPROCESSORS_ONLN);
        for(long i = 0;i<numCPU;i++) {
        }
    }

    std::vector<T> map(std::vector<boost::function<T ()> > funcs) {
        for(long i = 0;i<funcs.size();i++) {
        }
    }
private:
    boost::barrier _barrier;
    std::vector<boost::thread>;
};

