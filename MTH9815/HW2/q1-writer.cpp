#include <iostream>
#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>

using namespace boost::interprocess;

int main(int, char**) {

    try{
        // Create a shared memory object
        shared_memory_object shm(create_only, "SharedMemo", read_write); 

        // Set the size of the shared memory
        shm.truncate(sizeof(int));

        // Map the shared memory into this process's address space
        mapped_region region(shm, read_write);

        // Get a pointer to the shared memory
        int* data = static_cast<int*>(region.get_address());

        // Publish an int value
        *data = 42;

        std::cout << "Integer written to shared memory. Published value: " << *data << std::endl;
    } catch (interprocess_exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}

