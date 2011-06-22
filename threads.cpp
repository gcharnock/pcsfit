
#include "threads.hpp"
#include <unistd.h>

using namespace boost;
using namespace std;

class TestThread {
    TestThread(long _id) id(id) {
    }
    operator()() {
        cout << id << ":doing work" << endl;
        sleep(1.0);
        cout << id << ":work done" << endl;
    }
    long id;
};


Multithreader::Multithreader() {
}

Multithreader::exec() {
    _barrier.wait();
    //Threads do work here
    _barrier.wait();
    //Results
}

