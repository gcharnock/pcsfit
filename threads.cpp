
#include "threads.hpp"
#include <unistd.h>
#include <iostream>

using namespace boost;
using namespace std;

class TestThread {
    TestThread(long _id)
		: id(id) {
    }
    void operator()() {
        cout << id << ":doing work" << endl;
        sleep(1.0);
        cout << id << ":work done" << endl;
    }
    long id;
};


