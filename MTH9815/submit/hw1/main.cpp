//
//  main.cpp
//  hw1
//
//  Created by Yicheng.Sun on 14/11/2023.
//

#include <iostream>
#include "LinkedList.h"
#include "DoublyLinkedList.h"
#include "HashTable.h"

// test function of LinkedList
void test_linkedlist() {
    LinkedList<int> myList;

    // Add elements to the list
    myList.Add(10);
    myList.Add(20);
    myList.Add(30);

    // Print the size of the list
    std::cout << "Size: " << myList.Size() << std::endl;

    // Test Iterator
    ListIterator<int> iterator = myList.Iterator();
    std::cout << "List elements: ";
    while (iterator.HasNext()) {
        std::cout << iterator.Next() << " ";
    }
    std::cout << std::endl;

    // Test insert method
    myList.Insert(15, 1);
    std::cout << "Size after insertion: " << myList.Size() << std::endl;
    std::cout << "List after insertion: ";
    iterator = myList.Iterator();
    while (iterator.HasNext()) {
        std::cout << iterator.Next() << " ";
    }
    std::cout << std::endl;
    
    // Test get method
    std::cout << "Element at index 1: " << myList.Get(1) << std::endl;
    
    // Test indexof method
    int index = myList.IndexOf(20);
    std::cout << "Index of 20: " << index << std::endl;

    // Test remove method
    int removedValue = myList.Remove(2);
    std::cout << "Removed element at index 2: " << removedValue << std::endl;

    std::cout << "Size after removal: " << myList.Size() << std::endl;
    std::cout << "List after removal: ";
    iterator = myList.Iterator();
    while (iterator.HasNext()) {
        std::cout << iterator.Next() << " ";
    }
    std::cout << std::endl;
}

void test_doublylinkedlist() {
    DoublyLinkedList<int> myDoublyList;

    // Add elements to the list
    myDoublyList.Add(10);
    myDoublyList.Add(20);
    myDoublyList.Add(30);

    // Print the size of the list
    std::cout << "Size: " << myDoublyList.Size() << std::endl;

    // Test Iterator (forward)
    ListIterator<int> iterator = myDoublyList.Iterator();
    std::cout << "DList elements: ";
    while (iterator.HasNext()) {
        std::cout << iterator.Next() << " ";
    }
    std::cout << std::endl;

    // Test other methods
    myDoublyList.Insert(15, 1);
    
    std::cout << "Size after insertion: " << myDoublyList.Size() << std::endl;

    std::cout << "DList after insertion: ";
    iterator = myDoublyList.Iterator();
    while (iterator.HasNext()) {
        std::cout << iterator.Next() << " ";
    }
    std::cout << std::endl;

    std::cout << "Element at index 1: " << myDoublyList.Get(1) << std::endl;

    int index = myDoublyList.IndexOf(20);
    std::cout << "Index of 20: " << index << std::endl;

    int removedValue = myDoublyList.Remove(2);
    std::cout << "Removed element at index 2: " << removedValue << std::endl;

    // Print the size of the list after removal
    std::cout << "Size after removal: " << myDoublyList.Size() << std::endl;
    
    std::cout << "DList after removal: ";
    iterator = myDoublyList.Iterator();
    while (iterator.HasNext()) {
        std::cout << iterator.Next() << " ";
    }
    std::cout << std::endl;
}

void test_hashtable() {
    IntHasher intHasher;
    IntEqualityPredicate intEqualityPredicate;
    Hashtable<int, std::string> intHashtable(5, intHasher, intEqualityPredicate);

    // Insert some key-value pairs
    intHashtable.insert(1, "One");
    intHashtable.insert(2, "Two");
    intHashtable.insert(3, "Three");

    // Find and print values
    std::cout << "Try to find key: " << std::endl;
    std::string result;
    try {
        intHashtable.find(2, result);
        std::cout << "Value found for key = 2: " << result << std::endl;
    } catch (const std::runtime_error& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    
    // Remove a key
    intHashtable.remove(2);
    
    // Try to find the removed key
    std::cout << "After removal: " << std::endl;
    try {
        intHashtable.find(2, result);
        std::cout << "Value found for key = 2: " << result << std::endl;
    } catch (const std::runtime_error& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    
    // Try to insert a duplicate key (this should throw an error)
    std::cout << "Try to insert duplicate key. " << std::endl;
    try {
        intHashtable.insert(1, "DuplicateOne");
    } catch (const std::invalid_argument& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }

    // Clear the hashtable
    intHashtable.clear();
}

int main() {
    std::cout << "============Exercise 1============" << std::endl;
    test_linkedlist();
    
    std::cout << "============Exercise 2============" << std::endl;
    test_doublylinkedlist();
    
    std::cout << "============Exercise 3============" << std::endl;
    test_hashtable();

    return 0;
}

