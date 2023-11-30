/**
 * Definition of our Service Oriented Architecture (SOA) Service base class
 */

#ifndef SOA_HPP
#define SOA_HPP

/**
 * Definition of a generic base class Service.
 * Uses key generic type K and value generic type V.
 */
template<typename K, typename V>
class Service
{
public:
  // 
  virtual V& GetData(K key) = 0;
};

#endif
