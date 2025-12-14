#pragma once

template <typename T>
struct TernaryVector {
    int n;
    T values;
    T signs;
    bool valid;

    TernaryVector();
    TernaryVector(int n);
    TernaryVector(int n, int index);
    TernaryVector(int n, int *values);

    void set(int index, int value);
    void inverse();

    int compare(const TernaryVector &vector) const;

    bool operator==(const TernaryVector &vector) const;
    bool operator!=(const TernaryVector &vector) const;
    int operator[](int index) const;

    TernaryVector operator+(const TernaryVector &vector) const;
    TernaryVector operator-(const TernaryVector &vector) const;
    TernaryVector operator-() const;

    TernaryVector& operator+=(const TernaryVector &vector);
    TernaryVector& operator-=(const TernaryVector &vector);

    bool limit(bool checkFirstNonZero) const;
    bool limitSum(const TernaryVector &vector, bool checkFirstNonZero) const;
    bool limitSub(const TernaryVector &vector, bool checkFirstNonZero) const;

    bool positiveFirstNonZero() const;
    bool positiveFirstNonZeroSub(const TernaryVector &vector) const;

    int nonZeroCount() const;

    operator bool() const;

    template <typename T1>
    friend std::ostream& operator<<(std::ostream &os, const TernaryVector<T1> &vector);

    template <typename T1>
    friend std::istream& operator>>(std::ostream &is, TernaryVector<T1> &vector);
};

template <typename T>
TernaryVector<T>::TernaryVector() {
    n = 0;
    values = 0;
    signs = 0;
    valid = true;
}

template <typename T>
TernaryVector<T>::TernaryVector(int n) {
    this->n = n;
    this->values = 0;
    this->signs = 0;
    this->valid = true;
}

template <typename T>
TernaryVector<T>::TernaryVector(int n, int index) {
    this->n = n;
    this->values = T(1) << index;
    this->signs = 0;
    this->valid = true;
}

template <typename T>
TernaryVector<T>::TernaryVector(int n, int *values) {
    this->n = n;
    this->values = 0;
    this->signs = 0;
    this->valid = true;

    for (int i = 0; i < n; i++)
        set(i, values[i]);
}

template <typename T>
void TernaryVector<T>::set(int index, int value) {
    T mask = T(1) << index;

    if (value == 0) {
        values &= ~mask;
        signs &= ~mask;
    }
    else if (value == 1) {
        values |= mask;
        signs &= ~mask;
    }
    else if (value == -1) {
        values |= mask;
        signs |= mask;
    }
    else {
        valid = false;
        printf("invalid set (%d, %d)\n", index, value);
    }
}

template <typename T>
void TernaryVector<T>::inverse() {
    signs = (~signs) & values;
}

template <typename T>
int TernaryVector<T>::compare(const TernaryVector<T> &vector) const {
    if (values != vector.values)
        return 0;

    if (signs == vector.signs)
        return 1;

    if (signs == ((~vector.signs) & vector.values))
        return -1;

    return 0;
}

template <typename T>
bool TernaryVector<T>::operator==(const TernaryVector<T> &vector) const {
    return values == vector.values && signs == vector.signs;
}

template <typename T>
bool TernaryVector<T>::operator!=(const TernaryVector<T> &vector) const {
    return values != vector.values || signs != vector.signs;
}

template <typename T>
int TernaryVector<T>::operator[](int index) const {
    int value = int((values >> index) & 1);

    if ((signs >> index) & 1)
        value = -value;

    return value;
}

template <typename T>
TernaryVector<T> TernaryVector<T>::operator+(const TernaryVector<T> &vector) const {
    TernaryVector<T> result(n);

    result.values = values ^ vector.values;
    result.signs = ((signs & values) | (vector.signs & vector.values)) & result.values;
    result.valid = !(values & vector.values & ~(signs ^ vector.signs));
    return result;
}

template <typename T>
TernaryVector<T> TernaryVector<T>::operator-(const TernaryVector<T> &vector) const {
    TernaryVector<T> result(n);

    result.values = values ^ vector.values;
    result.signs = ((signs & values) | (~vector.signs & vector.values)) & result.values;
    result.valid = !(values & vector.values & (signs ^ vector.signs));
    return result;
}

template <typename T>
TernaryVector<T> TernaryVector<T>::operator-() const {
    TernaryVector<T> result(n);
    result.values = values;
    result.signs = (~signs) & values;
    result.valid = valid;
    return result;
}

template <typename T>
TernaryVector<T>& TernaryVector<T>::operator+=(const TernaryVector<T> &vector) {
    T sv1 = signs & values;
    T sv2 = vector.signs & vector.values;

    valid = !((values & vector.values & ~(signs ^ vector.signs)));
    values ^= vector.values;
    signs = (sv1 | sv2) & values;
    return *this;
}

template <typename T>
TernaryVector<T>& TernaryVector<T>::operator-=(const TernaryVector<T> &vector) {
    T sv1 = signs & values;
    T sv2 = ~vector.signs & vector.values;

    valid = !(values & vector.values & (signs ^ vector.signs));
    values ^= vector.values;
    signs = (sv1 | sv2) & values;
    return *this;
}

template <typename T>
TernaryVector<T>::operator bool() const {
    return values != 0;
}

template <typename T>
bool TernaryVector<T>::limit(bool checkFirstNonZero) const {
    if (!valid)
        return false;

    if (checkFirstNonZero)
        return values == 0 || (values & ~(values & (values - 1)) & ~signs);

    return true;
}

template <typename T>
bool TernaryVector<T>::limitSum(const TernaryVector<T> &vector, bool checkFirstNonZero) const {
    bool invalid = values & vector.values & ~(signs ^ vector.signs);

    if (invalid)
        return false;

    if (checkFirstNonZero) {
        T sumValues = values ^ vector.values;
        T sumSigns = ((signs & values) | (vector.signs & vector.values)) & sumValues;

        return sumValues == 0 || (sumValues & ~(sumValues & (sumValues - 1)) & ~sumSigns);
    }

    return true;
}

template <typename T>
bool TernaryVector<T>::limitSub(const TernaryVector<T> &vector, bool checkFirstNonZero) const {
    bool invalid = (values & vector.values & (signs ^ vector.signs));

    if (invalid)
        return false;

    if (checkFirstNonZero) {
        T subValues = values ^ vector.values;
        T subSigns = ((signs & values) | (~vector.signs & vector.values)) & subValues;

        return subValues == 0 || (subValues & ~(subValues & (subValues - 1)) & ~subSigns);
    }

    return true;
}

template <typename T>
bool TernaryVector<T>::positiveFirstNonZeroSub(const TernaryVector<T> &vector) const {
    T subValues = values ^ vector.values;
    T subSigns = ((signs & values) | (~vector.signs & vector.values)) & subValues;
    return subValues == 0 || (subValues & ~(subValues & (subValues - 1)) & ~subSigns);
}

template <typename T>
bool TernaryVector<T>::positiveFirstNonZero() const {
    return values == 0 || (values & ~(values & (values - 1)) & ~signs);
}

template <typename T>
int TernaryVector<T>::nonZeroCount() const {
    return __builtin_popcountll(values);
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const TernaryVector<T> &vector) {
    for (int i = 0; i < vector.n; i++) {
        if (i > 0)
            os << ", ";

        os << vector[i];
    }

    return os;
}

template <typename T>
std::istream& operator>>(std::istream &is, TernaryVector<T> &vector) {
    int value;
    vector.values = 0;
    vector.signs = 0;
    vector.valid = true;

    for (int i = 0; i < vector.n; i++) {
        is >> value;
        vector.set(i, value);
    }

    return is;
}
