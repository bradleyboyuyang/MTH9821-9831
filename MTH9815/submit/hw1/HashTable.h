//
//  HashTable.h
//  hw1
//
//  Created by Yicheng.Sun on 14/11/2023.
//

#ifndef HashTable_h
#define HashTable_h

#include <vector>
#include <list>
#include <stdexcept>

// Hasher class (base class)
template <typename K>
class Hasher {
public:
    virtual size_t hash(const K& key) const = 0;
};

// EqualityPredicate class (base class)
template <typename K>
class EqualityPredicate {
public:
    virtual bool isEqual(const K& key1, const K& key2) const = 0;
};

// Hashtable class
template <typename K, typename V>
class Hashtable {
private:
    size_t size;
    std::vector<std::list<std::pair<K, V>>> table;
    const Hasher<K>& hasher;
    const EqualityPredicate<K>& equalityPredicate;

public:
    // Constructor of Hashtable
    Hashtable(size_t initialSize, const Hasher<K>& hasher, const EqualityPredicate<K>& equalityPredicate)
        : size(initialSize), table(initialSize), hasher(hasher), equalityPredicate(equalityPredicate) {}
    
    // Insert a key-value pair into the hashtable
    void insert(const K& key, const V& value) {
        size_t index = hasher.hash(key) % size;
        
        // Check if the key already exists in the hash table
        for (const auto& entry : table[index]) {
            if (equalityPredicate.isEqual(entry.first, key)) {
                // Throw an exception for duplicate key
                throw std::invalid_argument("Duplicate key");
            }
        }
        
        // Insert the key-value pair into the hash table
        table[index].emplace_back(key, value);
    }
    
    // Find the value associated with the given key in the hashtable
    V& find(const K& key, V& value) const {
        size_t index = hasher.hash(key) % size;
        
        // Search for the key in the hash table
        for (const auto& entry : table[index]) {
            if (equalityPredicate.isEqual(entry.first, key)) {
                value = entry.second;
                return value;
            }
        }
        
        // Throw an exception if the key is not found
        throw std::runtime_error("Non-existing Key");
    }

    // Remove entries with the specified key from the hash table
    void remove(const K& key) {
        size_t index = hasher.hash(key) % size;
        auto& entryList = table[index];

        // Use std::remove_if to remove entries with the specified key
        auto it = std::remove_if(entryList.begin(), entryList.end(),
            [&](const auto& entry) { return equalityPredicate.isEqual(entry.first, key); });

        // Erase the removed entries from the list
        entryList.erase(it, entryList.end());
    }
    
    // Clear all entries from the hash table
    void clear() {
        for (auto& entryList : table) {
            entryList.clear();
        }
    }
};

// Example concrete Hasher class for integers
class IntHasher : public Hasher<int> {
public:
    size_t hash(const int& key) const override {
        return static_cast<size_t>(key);
    }
};

// Example concrete EqualityPredicate class for integers
class IntEqualityPredicate : public EqualityPredicate<int> {
public:
    bool isEqual(const int& key1, const int& key2) const override {
        return key1 == key2;
    }
};

#endif /* HashTable_h */
