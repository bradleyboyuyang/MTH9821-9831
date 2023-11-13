//
//  ListIterator.cpp
//  LinkedList
//
//  Created by Yicheng.Sun on 11/11/2023.
//

#include "ListIterator.hpp"
#include <stdexcept>

template <typename T>
ListIterator<T>::ListIterator(Node<T>* start) : current(start) {}

template <typename T>
bool ListIterator<T>::HasNext() {
    return current != nullptr;
}

template <typename T>
T& ListIterator<T>::Next() {
    // Throw an exception for invalid access
    if (current == nullptr) {
        throw std::out_of_range("Trying to access an element beyond the end of the list");
    }

    T& value = current->data;
    current = current->next;
    return value;
}
