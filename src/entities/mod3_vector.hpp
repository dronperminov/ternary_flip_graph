#pragma once

template <typename T>
struct Mod3Vector {
    int n;
    T low;
    T high;

    Mod3Vector();
    Mod3Vector(int n);
    Mod3Vector(int n, int index);
    Mod3Vector(int n, int *values);

    void set(int index, int value);
    void inverse();

    int compare(const Mod3Vector &vector) const;
    bool operator==(const Mod3Vector &vector) const;
    bool operator!=(const Mod3Vector &vector) const;
    int operator[](int index) const;

    Mod3Vector operator+(const Mod3Vector &vector) const;
    Mod3Vector operator-(const Mod3Vector &vector) const;
    Mod3Vector operator-() const;

    Mod3Vector& operator+=(const Mod3Vector &vector);
    Mod3Vector& operator-=(const Mod3Vector &vector);

    int nonZeroCount() const;

    operator bool() const;

    template <typename T1>
    friend std::ostream& operator<<(std::ostream &os, const Mod3Vector<T1> &vector);

    template <typename T1>
    friend std::istream& operator>>(std::ostream &is, Mod3Vector<T1> &vector);
};

template <typename T>
Mod3Vector<T>::Mod3Vector() {
    n = 0;
    low = 0;
    high = 0;
}

template <typename T>
Mod3Vector<T>::Mod3Vector(int n) {
    this->n = n;
    this->low = 0;
    this->high = 0;
}

template <typename T>
Mod3Vector<T>::Mod3Vector(int n, int index) {
    this->n = n;
    this->low = T(1) << index;
    this->high = 0;
}

template <typename T>
Mod3Vector<T>::Mod3Vector(int n, int *values) {
    this->n = n;
    this->low = 0;
    this->high = 0;

    for (int i = 0; i < n; i++)
        set(i, values[i]);
}

template <typename T>
void Mod3Vector<T>::set(int index, int value) {
    T mask = T(1) << index;

    value %= 3;
    if (value < 0)
        value += 3;

    if (value & 1)
        low |= mask;
    else
        low &= ~mask;

    if (value & 2)
        high |= mask;
    else
        high &= ~mask;
}

template <typename T>
void Mod3Vector<T>::inverse() {
    std::swap(low, high);
}

template <typename T>
int Mod3Vector<T>::compare(const Mod3Vector &vector) const {
    if (low == vector.low && high == vector.high)
        return 1;

    if (low == vector.high && high == vector.low)
        return -1;

    return 0;
}

template <typename T>
bool Mod3Vector<T>::operator==(const Mod3Vector<T> &vector) const {
    return low == vector.low && high == vector.high;
}

template <typename T>
bool Mod3Vector<T>::operator!=(const Mod3Vector<T> &vector) const {
    return low != vector.low || high != vector.high;
}

template <typename T>
int Mod3Vector<T>::operator[](int index) const {
    T mask = T(1) << index;
    return ((low & mask) ? 1 : 0) + ((high & mask) ? 2 : 0);
}

template <typename T>
Mod3Vector<T> Mod3Vector<T>::operator+(const Mod3Vector<T> &vector) const {
    Mod3Vector<T> result(n);
    T mask = (low | vector.low) & (high | vector.high);

    result.low = (low ^ vector.low) ^ (high & vector.high) ^ mask;
    result.high = (high ^ vector.high) ^ (low & vector.low) ^ mask;
    return result;
}

template <typename T>
Mod3Vector<T> Mod3Vector<T>::operator-(const Mod3Vector<T> &vector) const {
    Mod3Vector<T> result(n);
    T mask = (low | vector.high) & (high | vector.low);

    result.low = (low ^ vector.high) ^ (high & vector.low) ^ mask;
    result.high = (high ^ vector.low) ^ (low & vector.high) ^ mask;
    return result;
}

template <typename T>
Mod3Vector<T> Mod3Vector<T>::operator-() const {
    Mod3Vector<T> result(n);
    result.low = high;
    result.high = low;
    return result;
}

template <typename T>
Mod3Vector<T>& Mod3Vector<T>::operator+=(const Mod3Vector<T> &vector) {
    T mask = (low | vector.low) & (high | vector.high);
    T tmp = low;

    low = (low ^ vector.low) ^ (high & vector.high) ^ mask;
    high = (high ^ vector.high) ^ (tmp & vector.low) ^ mask;
    return *this;
}

template <typename T>
Mod3Vector<T>& Mod3Vector<T>::operator-=(const Mod3Vector<T> &vector) {
    T mask = (low | vector.high) & (high | vector.low);
    T tmp = low;

    low = (low ^ vector.high) ^ (high & vector.low) ^ mask;
    high = (high ^ vector.low) ^ (tmp & vector.high) ^ mask;
    return *this;
}

template <typename T>
Mod3Vector<T>::operator bool() const {
    return low | high;
}

template <typename T>
int Mod3Vector<T>::nonZeroCount() const {
    return __builtin_popcountll(low | high);
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const Mod3Vector<T> &vector) {
    for (int i = 0; i < vector.n; i++) {
        if (i > 0)
            os << ", ";

        os << vector[i];
    }

    return os;
}

template <typename T>
std::istream& operator>>(std::istream &is, Mod3Vector<T> &vector) {
    int value;
    vector.low = 0;
    vector.high = 0;

    for (int i = 0; i < vector.n; i++) {
        is >> value;
        vector.set(i, value);
    }

    return is;
}
