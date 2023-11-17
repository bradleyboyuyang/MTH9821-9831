//
//  LinkedList.h
//  hw1
//
//  Created by Yicheng.Sun on 14/11/2023.
//

#ifndef LinkedList_h
#define LinkedList_h

template <typename T>
class Node {
public:
    T data;
    Node<T>* next;

    Node(const T& value) : data(value), next(nullptr) {}
};

// forward declaration
template <typename T>
class ListIterator;

template <typename T>
class LinkedList {
protected:
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

};

template <typename T>
class ListIterator {
protected:
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

#endif /* LinkedList_h */
