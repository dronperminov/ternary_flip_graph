#pragma once

const int MIN_VALUE = -3;
const int MAX_VALUE = 3;

template <typename T>
struct IntegerVector {
    int n;
    std::vector<T> values;
    bool valid;

    IntegerVector();
    IntegerVector(int n);
    IntegerVector(int n, int index);
    IntegerVector(int n, int *values);

    void set(int index, int value);
    void inverse();

    int compare(const IntegerVector &vector) const;

    bool operator==(const IntegerVector &vector) const;
    bool operator!=(const IntegerVector &vector) const;
    int operator[](int index) const;

    IntegerVector operator+(const IntegerVector &vector) const;
    IntegerVector operator-(const IntegerVector &vector) const;
    IntegerVector operator-() const;

    IntegerVector& operator+=(const IntegerVector &vector);
    IntegerVector& operator-=(const IntegerVector &vector);

    bool limit(bool checkFirstNonZero) const;
    bool limitSum(const IntegerVector &vector, bool checkFirstNonZero) const;
    bool limitSub(const IntegerVector &vector, bool checkFirstNonZero) const;

    bool positiveFirstNonZero() const;
    bool positiveFirstNonZeroSub(const IntegerVector &vector) const;

    int nonZeroCount() const;

    operator bool() const;

    template <typename T1>
    friend std::ostream& operator<<(std::ostream &os, const IntegerVector<T1> &vector);

    template <typename T1>
    friend std::istream& operator>>(std::ostream &is, IntegerVector<T1> &vector);
};

template <typename T>
IntegerVector<T>::IntegerVector() {
    n = 0;
    valid = true;
}

template <typename T>
IntegerVector<T>::IntegerVector(int n) {
    this->n = n;
    this->values = std::vector<T>(n, 0);
    this->valid = true;
}

template <typename T>
IntegerVector<T>::IntegerVector(int n, int index) {
    this->n = n;
    this->values = std::vector<T>(n, 0);
    this->values[index] = 1;
    this->valid = true;
}

template <typename T>
IntegerVector<T>::IntegerVector(int n, int *values) {
    this->n = n;
    this->values = std::vector<T>(n, 0);
    this->valid = true;

    for (int i = 0; i < n; i++)
        set(i, values[i]);
}

template <typename T>
void IntegerVector<T>::set(int index, int value) {

    if (value >= MIN_VALUE && value <= MAX_VALUE) {
        this->values[index] = value;
    }
    else {
        valid = false;
        std::cout << "invalid set (" << index << ", " << value << ")" << std::endl;
    }
}

template <typename T>
void IntegerVector<T>::inverse() {
    for (int i = 0; i < n; i++)
        values[i] = -values[i];
}

template <typename T>
int IntegerVector<T>::compare(const IntegerVector<T> &vector) const {
    bool pos = true;
    bool neg = true;

    for (int i = 0; i < n && (pos || neg); i++) {
        pos &= values[i] == vector.values[i];
        neg &= values[i] == -vector.values[i];
    }

    if (pos)
        return 1;

    if (neg)
        return -1;

    return 0;
}

template <typename T>
bool IntegerVector<T>::operator==(const IntegerVector<T> &vector) const {
    for (int i = 0; i < n; i++)
        if (values[i] != vector.values[i])
            return false;

    return true;
}

template <typename T>
bool IntegerVector<T>::operator!=(const IntegerVector<T> &vector) const {
    for (int i = 0; i < n; i++)
        if (values[i] != vector.values[i])
            return true;

    return false;
}

template <typename T>
int IntegerVector<T>::operator[](int index) const {
    return values[index];
}

template <typename T>
IntegerVector<T> IntegerVector<T>::operator+(const IntegerVector<T> &vector) const {
    IntegerVector<T> result(n);

    for (int i = 0; i < n; i++) {
        result.values[i] = values[i] + vector.values[i];

        if (result.values[i] < MIN_VALUE || result.values[i] > MAX_VALUE)
            result.valid = false;
    }

    return result;
}

template <typename T>
IntegerVector<T> IntegerVector<T>::operator-(const IntegerVector<T> &vector) const {
    IntegerVector<T> result(n);

    for (int i = 0; i < n; i++) {
        result.values[i] = values[i] - vector.values[i];

        if (result.values[i] < MIN_VALUE || result.values[i] > MAX_VALUE)
            result.valid = false;
    }

    return result;
}

template <typename T>
IntegerVector<T> IntegerVector<T>::operator-() const {
    IntegerVector<T> result(n);

    for (int i = 0; i < n; i++) {
        result.values[i] = -values[i];

        if (result.values[i] < MIN_VALUE || result.values[i] > MAX_VALUE)
            result.valid = false;
    }

    return result;
}

template <typename T>
IntegerVector<T>& IntegerVector<T>::operator+=(const IntegerVector<T> &vector) {
    for (int i = 0; i < n; i++) {
        values[i] += vector.values[i];

        if (values[i] < MIN_VALUE || values[i] > MAX_VALUE)
            valid = false;
    }

    return *this;
}

template <typename T>
IntegerVector<T>& IntegerVector<T>::operator-=(const IntegerVector<T> &vector) {
    for (int i = 0; i < n; i++) {
        values[i] -= vector.values[i];

        if (values[i] < MIN_VALUE || values[i] > MAX_VALUE)
            valid = false;
    }

    return *this;
}

template <typename T>
IntegerVector<T>::operator bool() const {
    for (int i = 0; i < n; i++)
        if (values[i])
            return true;

    return false;
}

template <typename T>
bool IntegerVector<T>::limit(bool checkFirstNonZero) const {
    if (!valid)
        return false;

    if (checkFirstNonZero)
        for (int i = 0; i < n; i++)
            if (values[i] != 0)
                return values[i] > 0;

    return true;
}

template <typename T>
bool IntegerVector<T>::limitSum(const IntegerVector<T> &vector, bool checkFirstNonZero) const {
    T firstNonZero = 0;

    for (int i = 0; i < n; i++) {
        T value = values[i] + vector.values[i];

        if (value < MIN_VALUE || value > MAX_VALUE)
            return false;

        if (value != 0 && !firstNonZero)
            firstNonZero = value;
    }

    return !checkFirstNonZero || firstNonZero > 0;
}

template <typename T>
bool IntegerVector<T>::limitSub(const IntegerVector<T> &vector, bool checkFirstNonZero) const {
    T firstNonZero = 0;

    for (int i = 0; i < n; i++) {
        T value = values[i] - vector.values[i];

        if (value < MIN_VALUE || value > MAX_VALUE)
            return false;

        if (value != 0 && !firstNonZero)
            firstNonZero = value;
    }

    return !checkFirstNonZero || firstNonZero > 0;
}

template <typename T>
bool IntegerVector<T>::positiveFirstNonZeroSub(const IntegerVector<T> &vector) const {
    for (int i = 0; i < n; i++) {
        T value = values[i] - vector.values[i];

        if (value != 0)
            return value > 0;
    }

    return true;
}

template <typename T>
bool IntegerVector<T>::positiveFirstNonZero() const {
    for (int i = 0; i < n; i++)
        if (values[i] != 0)
            return values[i] > 0;

    return true;
}

template <typename T>
int IntegerVector<T>::nonZeroCount() const {
    int count = 0;

    for (int i = 0; i < n; i++)
        if (values[i] != 0)
            count++;

    return count;
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const IntegerVector<T> &vector) {
    for (int i = 0; i < vector.n; i++) {
        if (i > 0)
            os << ", ";

        os << vector.values[i];
    }

    return os;
}

template <typename T>
std::istream& operator>>(std::istream &is, IntegerVector<T> &vector) {
    vector.valid = true;

    for (int i = 0; i < vector.n; i++) {
        int value;
        is >> value;
        vector.set(i, value);
    }

    return is;
}
