//
//  ListIterator.hpp
//  LinkedList
//
//  Created by Yicheng.Sun on 11/11/2023.
//

#ifndef ListIterator_hpp
#define ListIterator_hpp

#include <stdio.h>
#include <stdexcept>
#include "Node.hpp"

// Forward declaration
template <typename T>
class LinkedList;


template <typename T>
// Iterator class for LinkedList
class ListIterator {
private:
    Node<T>* current;

public:
    // constructor to initialize an iterator with a starting node
    ListIterator(Node<T>* startNode) : current(startNode) {};

    // check if there is another element to return in this iterator
    bool HasNext() {
        return current != nullptr;
    };

    // return the next element in this iterator
    T& Next() {
        // Throw an exception for invalid access
        if (current == nullptr) {
            throw std::out_of_range("Trying to access an element beyond the end of the list");
        }

        T& value = current->data;
        current = current->next;
        return value;
    };
};

#endif /* ListIterator_hpp */
