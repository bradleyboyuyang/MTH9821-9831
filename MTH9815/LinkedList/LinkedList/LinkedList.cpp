//
//  LinkedList.cpp
//  LinkedList
//
//  Created by Yicheng.Sun on 11/11/2023.
//

#include "LinkedList.hpp"
#include "iostream"

template <typename T>
LinkedList<T>::LinkedList() : head(nullptr), size(0) {}

//template <typename T>
//void LinkedList<T>::Add(const T& value) {
//    Node<T>* newNode = new Node<T>(value);
//    
//    // check if the LinkedList is empty
//    if (head == nullptr) {
//        head = newNode;
//    } else {
//        Node<T>* current = head;
//        while (current->next != nullptr) {
//            current = current->next;
//        }
//        current->next = newNode;
//    }
//    
//    // update the LinkedList size
//    size++;
//}
//
//template <typename T>
//void LinkedList<T>::Insert(const T& value, int index) {
//    // check the validity of index
//    if (index < 0 || index > size) {
//        std::cerr << "Index out of range." << std::endl;
//        return;
//    }
//
//    Node<T>* newNode = new Node<T>(value);
//
//    if (index == 0) {
//        newNode->next = head;
//        head = newNode;
//    } else {
//        Node<T>* current = head;
//        for (int i = 0; i < index - 1; i++) {
//            current = current->next;
//        }
//        newNode->next = current->next;
//        current->next = newNode;
//    }
//
//    // update the LinkedList size
//    size++;
//}
//
//template <typename T>
//T& LinkedList<T>::Get(int index) {
//    // check the validity of index
//    if (index < 0 || index >= size) {
//        throw std::out_of_range("Index out of range.");
//    }
//
//    Node<T>* current = head;
//    for (int i = 0; i < index; i++) {
//        current = current->next;
//    }
//
//    return current->data;
//}
//
//template <typename T>
//int LinkedList<T>::IndexOf(const T& value) {
//    int index = 0;
//    Node<T>* current = head;
//    while (current != nullptr) {
//        if (current->data == value) {
//            return index;
//        }
//        current = current->next;
//        index++;
//    }
//
//    return -1;
//}
//
//template <typename T>
//T& LinkedList<T>::Remove(int index) {
//    // check the validity of index
//    if (index < 0 || index >= size) {
//        throw std::out_of_range("Index out of range.");
//    }
//
//    Node<T>* removedNode;
//    if (index == 0) {
//        removedNode = head;
//        head = head->next;
//    } else {
//        Node<T>* current = head;
//        for (int i = 0; i < index - 1; i++) {
//            current = current->next;
//        }
//        removedNode = current->next;
//        current->next = current->next->next;
//    }
//
//    T& removedData = removedNode->data;
//    delete removedNode;
//    
//    // update the LinkedList size
//    size--;
//
//    return removedData;
//}
//
//template <typename T>
//int LinkedList<T>::Size() {
//    return size;
//}
//
//template <typename T>
//ListIterator<T> LinkedList<T>::Iterator() {
//    return ListIterator<T>(head);
//}
//
//template <typename T>
//LinkedList<T>::~LinkedList() {
//    Node<T>* current = head;
//    Node<T>* next;
//    while (current != nullptr) {
//        next = current->next;
//        delete current;
//        current = next;
//    }
//}
