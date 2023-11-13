//
//  Node.hpp
//  LinkedList
//
//  Created by Yicheng.Sun on 11/11/2023.
//

#ifndef Node_hpp
#define Node_hpp

#include <stdio.h>

template <typename T>

// Node class for LinkedList
class Node {
public:
    T data;
    Node<T>* next;
    
    // constructor to initialize a Node with given value
    Node(const T& value) : data(value), next(nullptr) {};
};

#endif /* Node_hpp */
