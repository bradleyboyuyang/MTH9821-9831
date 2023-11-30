#pragma once
#include <vector>
#include <list>
#include <stdexcept>
#include <algorithm>

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
    Hashtable(size_t initialSize, const Hasher<K>& hasher, const EqualityPredicate<K>& equalityPredicate)
        : size(initialSize), table(initialSize), hasher(hasher), equalityPredicate(equalityPredicate) {}

    void insert(const K& key, const V& value) {
        size_t index = hasher.hash(key) % size;
        for (const auto& entry : table[index]) {
            if (equalityPredicate.isEqual(entry.first, key)) {
                throw std::invalid_argument("Duplicate key");
            }
        }
        table[index].emplace_back(key, value);
    }

    V& find(const K& key, V& value) const {
        size_t index = hasher.hash(key) % size;
        for (const auto& entry : table[index]) {
            if (equalityPredicate.isEqual(entry.first, key)) {
                value = entry.second;
                return value;
            }
        }
        throw std::runtime_error("Non-existing Key");
    }

    void remove(const K& key) {
        size_t index = hasher.hash(key) % size;
        auto& entryList = table[index];

        auto it = std::remove_if(entryList.begin(), entryList.end(),
            [&](const auto& entry) { return equalityPredicate.isEqual(entry.first, key); });

        entryList.erase(it, entryList.end());
    }

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