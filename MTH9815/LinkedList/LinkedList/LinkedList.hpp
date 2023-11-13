//
//  LinkedList.hpp
//  LinkedList
//
//  Created by Yicheng.Sun on 11/11/2023.
//

#ifndef LinkedList_hpp
#define LinkedList_hpp

#include <stdio.h>
#include "Node.hpp"
//#include "ListIterator.hpp"

template <typename T>

// LinkedList class with basic operations
class LinkedList {
private:
    Node<T>* head;
    int size;

public:
    // constructor to initialize an empty LinkedList
    LinkedList();

//    // add the specified element at the end of the list
//    void Add(const T& value);
//
//    // add the specified element at the specified index
//    void Insert(const T& value, int index);
//
//    // get the element at the specified index
//    T& Get(int index);
//
//    // retrieve the index of the specified element (-1 if it does not exist in the list)
//    int IndexOf(const T& value);
//
//    // remove the element at the specified index and return it
//    T& Remove(int index);
//
//    // return the size of the list
//    int Size();
//    
//    // return an iterator on this list
//    ListIterator<T> Iterator();
//
//    // destructor to clean up memory used by the list
//    ~LinkedList();
};

#endif /* LinkedList_hpp */
