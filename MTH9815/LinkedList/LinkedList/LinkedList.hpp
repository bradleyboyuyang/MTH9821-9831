//
//  LinkedList.hpp
//  LinkedList
//
//  Created by Yicheng.Sun on 11/11/2023.
//

#ifndef LinkedList_hpp
#define LinkedList_hpp

#include <stdio.h>
#include "iostream"
#include "Node.hpp"
#include "ListIterator.hpp"

template <typename T>

// LinkedList class with basic operations
class LinkedList {
private:
    Node<T>* head;
    int size;

public:
    // constructor to initialize an empty LinkedList
    LinkedList() : head(nullptr), size(0) {};

    // add the specified element at the end of the list
    void Add(const T& value) {
        Node<T>* newNode = new Node<T>(value);
        
        // check if the LinkedList is empty
        if (head == nullptr) {
            head = newNode;
        } else {
            Node<T>* current = head;
            while (current->next != nullptr) {
                current = current->next;
            }
            current->next = newNode;
        }
        
        // update the LinkedList size
        size++;
    };

    // add the specified element at the specified index
    void Insert(const T& value, int index) {
        // check the validity of index
        if (index < 0 || index > size) {
            std::cerr << "Index out of range." << std::endl;
            return;
        }

        Node<T>* newNode = new Node<T>(value);

        if (index == 0) {
            newNode->next = head;
            head = newNode;
        } else {
            Node<T>* current = head;
            for (int i = 0; i < index - 1; i++) {
                current = current->next;
            }
            newNode->next = current->next;
            current->next = newNode;
        }

        // update the LinkedList size
        size++;
    };

    // get the element at the specified index
    T& Get(int index) {
        // check the validity of index
        if (index < 0 || index >= size) {
            throw std::out_of_range("Index out of range.");
        }

        Node<T>* current = head;
        for (int i = 0; i < index; i++) {
            current = current->next;
        }

        return current->data;
    };

    // retrieve the index of the specified element (-1 if it does not exist in the list)
    int IndexOf(const T& value) {
        int index = 0;
        Node<T>* current = head;
        while (current != nullptr) {
            if (current->data == value) {
                return index;
            }
            current = current->next;
            index++;
        }

        return -1;
    };

    // remove the element at the specified index and return it
    T& Remove(int index){
        // check the validity of index
        if (index < 0 || index >= size) {
            throw std::out_of_range("Index out of range.");
        }

        Node<T>* removedNode;
        if (index == 0) {
            removedNode = head;
            head = head->next;
        } else {
            Node<T>* current = head;
            for (int i = 0; i < index - 1; i++) {
                current = current->next;
            }
            removedNode = current->next;
            current->next = current->next->next;
        }

        T& removedData = removedNode->data;
        
        // update the LinkedList size
        size--;

        return removedData;
    };

    // return the size of the list
    int Size() {
        return size;
    };
    
    // return an iterator on this list
    ListIterator<T> Iterator() {
        return ListIterator<T>(head);
    };

    // destructor to clean up memory used by the list
    ~LinkedList() {
        Node<T>* current = head;
        Node<T>* next;
        while (current != nullptr) {
            next = current->next;
            delete current;
            current = next;
        }
    };
};

#endif /* LinkedList_hpp */
