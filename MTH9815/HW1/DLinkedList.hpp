#pragma once
#include "LinkedList.hpp"

template <typename T>
class DNode : public Node<T> {
public:
    DNode<T>* prev;
    DNode(const T& value) : Node<T>(value), prev(nullptr) {}
};


template <typename T>
class DoublyLinkedList : public LinkedList<T> {
protected:
    DNode<T>* tail;

public:
    DoublyLinkedList() : LinkedList<T>(), tail(nullptr) {}

    void Add(const T& value) {
        DNode<T>* newNode = new DNode<T>(value);

        if (this->head == nullptr) {
            this->head = newNode;
        } else {
            tail->next = newNode;
            newNode->prev = tail;
        }
        
        tail = newNode;
        this->size++;
    };
    
    void Insert(const T& value, int index) {
        if (index < 0 || index > this->size) {
            std::cerr << "Index out of range." << std::endl;
            return;
        }

        DNode<T>* newNode = new DNode<T>(value);

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

        this->size++;
    };
    
    T& Remove(int index) {
        if (index < 0 || index >= this->size) {
            throw std::out_of_range("Index out of range.");
        }

        DNode<T>* removedNode;

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
        
        this->size--;

        return removedData;
    };
};

