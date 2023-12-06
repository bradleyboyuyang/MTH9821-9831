#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <iostream>

using namespace boost::interprocess;

int main() {
    try {
        // open shared memory object
        shared_memory_object shm(open_only, "SharedMemo", read_only);

        // map shared memory to address place
        mapped_region region(shm, read_only);

        // take the pointer and read the data
        const int* data = static_cast<int*>(region.get_address());
        std::cout << "Integer read from shared memory: " << *data << std::endl;

        // remove the shared memory
        shared_memory_object::remove("SharedMemo");
    } catch (interprocess_exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}
