//
//  DoublyLinkedList.h
//  hw1
//
//  Created by Yicheng.Sun on 14/11/2023.
//

#ifndef DoublyLinkedList_h
#define DoublyLinkedList_h

#include "LinkedList.h"


template <typename T>
class DNode : public Node<T> {
public:
    DNode<T>* prev;
    
    // Constructor of DNode
    DNode(const T& value) : Node<T>(value), prev(nullptr) {}
};


template <typename T>
class DoublyLinkedList : public LinkedList<T> {
protected:
    DNode<T>* tail;

public:
    // Constructor to initialize and empty DoublyLinkedList
    DoublyLinkedList() : LinkedList<T>(), tail(nullptr) {}
    
    // Add the specified element at the end of the list
    void Add(const T& value) {
        DNode<T>* newNode = new DNode<T>(value);
        
        // Check if the DoublyLinkedList is empty
        if (this->head == nullptr) {
            this->head = newNode;
        } else {
            tail->next = newNode;
            newNode->prev = tail;
        }
        
        tail = newNode;
        // update the DoublyLinkedList size
        this->size++;
    };
    
    // Add the specified element at the specified index
    void Insert(const T& value, int index) {
        // Check the validity of index
        if (index < 0 || index > this->size) {
            std::cerr << "Index out of range." << std::endl;
            return;
        }

        DNode<T>* newNode = new DNode<T>(value);

        // Check if the inserted DNode is the head
        if (index == 0) {
            newNode->next = this->head;
            if (this->head != nullptr) {
                static_cast<DNode<T>*>(this->head)->prev = newNode;
            }
            this->head = newNode;
        } else {
            DNode<T>* current = static_cast<DNode<T>*>(this->head);
            for (int i = 0; i < index - 1; i++) {
                current = static_cast<DNode<T>*>(current->next);
            }
            newNode->next = current->next;
            newNode->prev = current;
            if (current->next != nullptr) {
                static_cast<DNode<T>*>(current->next)->prev = newNode;
            }
            current->next = newNode;
        }

        // update the DoublyLinkedList size
        this->size++;
    };
    
    // Remove the element at the specified index and return it
    T& Remove(int index) {
        // Check the validity of index
        if (index < 0 || index >= this->size) {
            throw std::out_of_range("Index out of range.");
        }

        DNode<T>* removedNode;

        // Check if the removed DNode is the head
        if (index == 0) {
            removedNode = static_cast<DNode<T>*>(this->head);
            this->head = this->head->next;
            if (this->head != nullptr) {
                static_cast<DNode<T>*>(this->head)->prev = nullptr;
            }
        } else {
            DNode<T>* current = static_cast<DNode<T>*>(this->head);
            for (int i = 0; i < index - 1; i++) {
                current = static_cast<DNode<T>*>(current->next);
            }
            removedNode = static_cast<DNode<T>*>(current->next);
            current->next = current->next->next;
            if (current->next != nullptr) {
                static_cast<DNode<T>*>(current->next)->prev = current;
            }
        }

        T& removedData = removedNode->data;
        
        // update the DoublyLinkedList size
        this->size--;

        return removedData;
    };
};


#endif /* DoublyLinkedList_h */
