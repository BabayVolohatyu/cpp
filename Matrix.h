#ifndef MATRIX_H
#define MATRIX_H

#define MIN_LIMIT 1e-10

template <typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
class Matrix {
private:
    std::vector<std::vector<T> > matrix;

    static void approximate(std::vector<std::vector<T>>& matrix) {
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                if (std::abs(matrix[i][j]) <= MIN_LIMIT) { matrix[i][j] = 0; }
            }
        }
    }

public:
    Matrix(size_t rows = 1, size_t cols = 1, T initialValue = T()) {
        matrix.resize(rows, std::vector<T>(cols, initialValue));
    }

    Matrix(Matrix &&other) noexcept : matrix(std::move(other.matrix)) {
    }

    Matrix(const Matrix &other) : matrix(other.matrix) {
    }

    ~Matrix() = default;

    Matrix transpose() const {
        Matrix result(matrix[0].size(), matrix.size());
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                result[j][i] = matrix[i][j];
            }
        }
        return result;
    }

    Matrix inverse() {
        size_t n = matrix.size();
        if (n != matrix[0].size()) {
            throw std::invalid_argument("Matrix must be square to compute inverse.");
        }

        // Create an augmented matrix [A | I]
        Matrix augmented(n, 2 * n);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                augmented[i][j] = matrix[i][j];
                augmented[i][j + n] = (i == j) ? 1 : 0; // Identity matrix
            }
        }

        // Perform Gaussian elimination
        for (size_t i = 0; i < n; i++) {
            // Find the pivot row
            size_t pivot = i;
            for (size_t j = i + 1; j < n; j++) {
                if (std::abs(augmented[j][i]) > std::abs(augmented[pivot][i])) {
                    pivot = j;
                }
            }

            // Check for singular matrix
            if (augmented[pivot][i] == 0) {
                throw std::invalid_argument("Matrix is singular and cannot be inverted.");
            }

            // Normalize the pivot row
            T pivotValue = augmented[i][i];
            for (size_t j = 0; j < 2 * n; j++) {
                augmented[i][j] /= pivotValue;
            }

            // Eliminate other rows
            for (size_t j = 0; j < n; j++) {
                if (j != i) {
                    T factor = augmented[j][i];
                    for (size_t k = 0; k < 2 * n; k++) {
                        augmented[j][k] -= factor * augmented[i][k];
                    }
                }
            }
        }

        // Extract the inverse from the augmented matrix
        Matrix inverse(n, n);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                inverse[i][j] = augmented[i][j + n];
            }
        }
        approximate(augmented.matrix);
        return inverse;
    }

    void resize(size_t rows, size_t cols, T initialValue = T()) {
        matrix.resize(rows, std::vector<T>(cols, initialValue));
    }

    void clear() {
        matrix.clear();
    }

    std::vector<T> &operator[](size_t row) {
        return matrix[row];
    }

    Matrix &operator=(Matrix &&other) noexcept {
        if (this != &other) {
            matrix = std::move(other.matrix);
        }
        return *this;
    }

    Matrix &operator=(const Matrix &other) {
        if (this != &other) {
            matrix = other.matrix;
        }
        return *this;
    }

    Matrix operator*(const Matrix &other) {
        size_t numberOfRows = matrix.size(); // number of rows of the first matrix
        size_t numberOfColumns = other.size().second; // number of columns of the second matrix
        size_t n = matrix[0].size(); //columns of first matrix
        //check if multiplication is available
        if (n != other.size().first) {
            throw std::invalid_argument("Matrix dimensions don`t match");
        }

        Matrix result(numberOfRows, numberOfColumns); //result matrix

        for (size_t i = 0; i < numberOfRows; i++) {
            for (size_t j = 0; j < numberOfColumns; j++) {
                for (size_t k = 0; k < n; k++) {
                    result.matrix[i][j] += matrix[i][k] * other.matrix[k][j];
                }
            }
        }

        approximate(result.matrix);
        return result;
    }

    Matrix operator*(const T &scalar) {
        Matrix result(matrix.size(), matrix[0].size());
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                result.matrix[i][j] = matrix[i][j] * scalar;
            }
        }

        approximate(result.matrix);
        return result;
    }

    Matrix &operator*=(const Matrix &other) {
        // Check if multiplication is available
        if (matrix[0].size() != other.size().first) {
            throw std::invalid_argument("Matrix dimensions don't match");
        }

        // Create a temporary matrix to store the result
        Matrix result(matrix.size(), other.size().second);

        // Perform matrix multiplication
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < other.size().second; j++) {
                result.matrix[i][j] = 0; // Initialize the result element to 0
                for (size_t k = 0; k < matrix[0].size(); k++) {
                    result.matrix[i][j] += matrix[i][k] * other.matrix[k][j];
                }
            }
        }

        approximate(result.matrix);
        // Assign the result back to the current matrix
        matrix = std::move(result.matrix);

        return *this;
    }

    Matrix &operator*=(const T &scalar) {
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                matrix[i][j] *= scalar;
            }
        }

        approximate(matrix);
        return *this;
    }

    Matrix operator/(const T &scalar) {
        if (scalar == 0.0) {
            throw std::invalid_argument("Cannot divide by zero");
        }
        Matrix result(matrix.size(), matrix[0].size());
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                result.matrix[i][j] = matrix[i][j] / scalar;
            }
        }

        approximate(result.matrix);
        return result;
    }

    Matrix &operator/=(const T &scalar) {
        if (scalar == 0.0) {
            throw std::invalid_argument("Cannot divide by zero");
        }
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                matrix[i][j] /= scalar;
            }
        }

        approximate(matrix);
        return *this;
    }

    Matrix operator+(const Matrix &other) {
        // Check if dimensions match
        if (matrix.size() != other.size().first || matrix[0].size() != other.size().second) {
            throw std::invalid_argument("Matrix dimensions don't match");
        }

        Matrix result(matrix.size(), matrix[0].size()); // Result matrix

        // Perform element-wise addition
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                result.matrix[i][j] = matrix[i][j] + other.matrix[i][j];
            }
        }

        return result;
    }

    Matrix operator+(const T &scalar) {
        Matrix result(matrix.size(), matrix[0].size());
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                result.matrix[i][j] = matrix[i][j] + scalar;
            }
        }
        return result;
    }

    Matrix &operator+=(const Matrix &other) {
        // Check if dimensions match
        if (matrix.size() != other.size().first || matrix[0].size() != other.size().second) {
            throw std::invalid_argument("Matrix dimensions don't match");
        }

        // Perform element-wise addition and assign to the current matrix
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                matrix[i][j] += other.matrix[i][j];
            }
        }

        return *this;
    }

    Matrix &operator+=(const T &scalar) {
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                matrix[i][j] += scalar;
            }
        }
        return *this;
    }

    Matrix operator-(const Matrix &other) {
        // Check if dimensions match
        if (matrix.size() != other.size().first || matrix[0].size() != other.size().second) {
            throw std::invalid_argument("Matrix dimensions don't match");
        }

        Matrix result(matrix.size(), matrix[0].size()); // Result matrix

        // Perform element-wise subtraction
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                result.matrix[i][j] = matrix[i][j] - other.matrix[i][j];
            }
        }

        return result;
    }

    Matrix operator-(const T &scalar) {
        Matrix result(matrix.size(), matrix[0].size());
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                result.matrix[i][j] = matrix[i][j] - scalar;
            }
        }
        return result;
    }

    Matrix &operator-=(const Matrix &other) {
        // Check if dimensions match
        if (matrix.size() != other.size().first || matrix[0].size() != other.size().second) {
            throw std::invalid_argument("Matrix dimensions don't match");
        }

        // Perform element-wise subtraction and assign to the current matrix
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                matrix[i][j] -= other.matrix[i][j];
            }
        }

        return *this;
    }

    Matrix &operator-=(const T &scalar) {
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                matrix[i][j] -= scalar;
            }
        }
        return *this;
    }

    bool operator==(const Matrix &other) const {
        return matrix == other.matrix;
    }

    bool operator!=(const Matrix &other) const {
        return !(*this == other);
    }

    friend std::ostream &operator<<(std::ostream &os, const Matrix &obj) {
        for (size_t i = 0; i < obj.matrix.size(); i++) {
            os << "[ ";
            for (size_t j = 0; j < obj.matrix[0].size(); j++) {
                os << obj.matrix[i][j] << " ";
            }
            os << ']' << std::endl;
        }
        return os;
    }
    std::pair<size_t, size_t> size() const {
        return {matrix.size(), matrix[0].size()};
    }

    T get(size_t row, size_t col) const {
        if (row >= matrix.size() || col >= matrix[0].size()) {
            throw std::out_of_range("Matrix indices out of range");
        }
        return matrix[row][col];
    }

    void set(size_t row, size_t col, T value) {
        if (row >= matrix.size() || col >= matrix[0].size()) {
            throw std::out_of_range("Matrix indices out of range");
        }
        matrix[row][col] = value;
    }
};
#endif
