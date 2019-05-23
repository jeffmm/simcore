// Helper functions, such as enums, etc

#ifndef _SIMCORE_HELPERS_H_
#define _SIMCORE_HELPERS_H_

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <map>

// Taken from stack overflow
// http://codereview.stackexchange.com/questions/14309/conversion-between-enum-and-string-in-c-class-header
namespace enh {
    template<typename T>
    struct enumStrings {
        static char const* data[];
    };

    // Utility type
    template<typename T>
    struct enumRefHolder {
        T& enumVal;
        enumRefHolder(T& enumVal): enumVal(enumVal) {}
    };
    template<typename T>
    struct enumConstRefHolder {
        T const& enumVal;
        enumConstRefHolder(T const& enumVal): enumVal(enumVal) {}
    };

    // Actual work functions
    template<typename T>
    std::ostream& operator<<(std::ostream& str, enumConstRefHolder<T> const& data) {
        return str << enumStrings<T>::data[data.enumVal];
    }

    template<typename T>
    std::istream& operator>>(std::istream& str, enumRefHolder<T> const& data) {
        std::string value;
        str >> value;

        static auto begin   = std::begin(enumStrings<T>::data);
        static auto end     = std::end(enumStrings<T>::data);

        auto find = std::find(begin, end, value);
        if (find != end) {
            data.enumVal = static_cast<T>(std::distance(begin, find));
        }
        return str;
    }

    // Finally, the public interface
    template<typename T>
    enumConstRefHolder<T>   enumToString(T const& e) { return enumConstRefHolder<T>(e); }

    template<typename T>
    enumRefHolder<T>        enumFromString(T& e) { return enumRefHolder<T>(e); }
}

// Trying to do reflection in c++11, eww
namespace rfh {
    template <class T> void* constructor() { return (void*)new T(); }

    struct factory {
        typedef void*(*constructor_t)();
        typedef std::map<std::string, constructor_t> map_type;
        map_type m_classes;

        template<class T>
        void register_class(std::string const& n) {
            m_classes.insert(std::make_pair(n, &constructor<T>));
        }

        void* construct(std::string const& n) {
            map_type::iterator i = m_classes.find(n);
            if (i == m_classes.end()) return 0;
            return i->second();
        }

        void clear() {
          m_classes.clear();
        }
    };
}

// hash helpers
namespace hashh {
  template <class T>
  inline void hash_combine(std::size_t &seed, const T &v) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  struct pair_hash {
    template<typename S, typename T>
    inline std::size_t operator()(const std::pair<S, T> &v) const {
      std::size_t seed = 0;
      hash_combine(seed, v.first);
      hash_combine(seed, v.second);
      return seed;
    }
  };
}

// yaml helpers to find if something failed
namespace yamlh {
}

#endif
