//
//  ListIterator.hpp
//  LinkedList
//
//  Created by Yicheng.Sun on 11/11/2023.
//

#ifndef ListIterator_hpp
#define ListIterator_hpp

#include <stdio.h>
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
    ListIterator(Node<T>* startNode);

    // check if there is another element to return in this iterator
    bool HasNext();

    // return the next element in this iterator
    T& Next();
};

#endif /* ListIterator_hpp */
