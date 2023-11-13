//
//  main.cpp
//  LinkedList
//
//  Created by Yicheng.Sun on 11/11/2023.
//

#include <iostream>
#include "LinkedList.hpp"

int main() {
    // Create a LinkedList of integers
    LinkedList<int> myList;

    // Add elements to the list
    myList.Add(1);
    myList.Add(2);
    myList.Add(3);

    // Insert an element at index 1
    myList.Insert(4, 1);

    // Print the elements using an iterator
    ListIterator<int> iterator = myList.Iterator();
    std::cout << "List after insertion: ";
    while (iterator.HasNext()) {
        std::cout << iterator.Next() << " ";
    }
    std::cout << std::endl;

    // Get and print the size of the list
    std::cout << "Size of the list: " << myList.Size() << std::endl;

    // Get the element at index 2
    int elementAtIndex2 = myList.Get(2);
    std::cout << "Element at index 2: " << elementAtIndex2 << std::endl;

    // Get the index of an element
    int index = myList.IndexOf(3);
    std::cout << "Index of element 3: " << index << std::endl;
    
    // Remove the element at index 1
    int removedElement = myList.Remove(1);

    // Print the elements after removal
    iterator = myList.Iterator();
    std::cout << "List after removal: ";
    while (iterator.HasNext()) {
        std::cout << iterator.Next() << " ";
    }
    std::cout << std::endl;

    // Print the removed element
    std::cout << "Removed element: " << removedElement << std::endl;

    // Print the updated size of the list
    std::cout << "Size of the list after removal: " << myList.Size() << std::endl;

    return 0;
}
