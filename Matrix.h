#ifndef MATRIX_H
#define MATRIX_H

template<class T>
class Matrix {
private:
    std::vector<std::vector<T> > matrix;

    class RowProxy {
    public:
        std::vector<T> &row;
        size_t col;

        RowProxy(std::vector<T> &r, size_t c = 0) : row(r), col(c) {
        }

        RowProxy(const RowProxy &r) : row(r.row), col(r.col) {}

        RowProxy &operator=(const RowProxy &r) {
          if (this != &r) {
            row = r.row;
            col = r.col;
          }
          return *this;
        }
        T &operator[](size_t c) {
            col = c;
            return row[col];
        }

        // I didnt write +,-,*,/ because its unnecessary and useless
        // Compound assignment with scalar
        RowProxy &operator+=(const T &scalar) {
            row[col] += scalar;
            return *this;
        }

        RowProxy &operator-=(const T &scalar) {
            row[col] -= scalar;
            return *this;
        }

        RowProxy &operator*=(const T &scalar) {
            row[col] *= scalar;
            return *this;
        }

        RowProxy &operator/=(const T &scalar) {
            row[col] /= scalar;
            return *this;
        }
    };

public:
    Matrix(size_t rows = 1, size_t cols = 1) {
        matrix.resize(rows, std::vector<T>(cols, 0));
    }

    Matrix(Matrix &&other) noexcept : matrix(std::move(other.matrix)) {
    }

    Matrix transpose() const {
        std::vector<std::vector<T>> transposed(matrix[0].size(), std::vector<T>(matrix.size()));

        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                transposed[j][i] = matrix[i][j];
            }
        }

        Matrix result;
        result.matrix = transposed;
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
            if (augmented[i][i] == 0) {
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

        return inverse;
    }

    class ConstRowProxy {
    private:
        const std::vector<T> &row;

    public:
        ConstRowProxy(const std::vector<T> &r) : row(r) {
        }

        const T &operator[](size_t col) const {
            return row[col];
        }
    };

    RowProxy operator[](size_t row) {
        return RowProxy(matrix[row]);
    }

    ConstRowProxy operator[](size_t row) const {
        return ConstRowProxy(matrix[row]);
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

    Matrix operator*(const Matrix &other) const {
        size_t numberOfRows = matrix.size(); // number of rows of the first matrix
        size_t numberOfColumns = other.size().second; // number of columns of the second matrix
        size_t n = matrix[0].size(); //columns of first matrix
        //check if multiplication is available
        if (n != other.size().first) {
            throw std::invalid_argument("Matrix dimensions don`t match");
        }

        Matrix result(numberOfRows, numberOfColumns); //result matrix

        for (int i = 0; i < numberOfRows; i++) {
            for (int j = 0; j < numberOfColumns; j++) {
                for (int k = 0; k < n; k++) {
                    result.matrix[i][j] += matrix[i][k] * other.matrix[k][j];
                }
            }
        }

        return result;
    }

    Matrix operator*(const T &scalar) const {
        Matrix result(matrix.size(), matrix[0].size());
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                result.matrix[i][j] = matrix[i][j] * scalar;
            }
        }
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
        for (int i = 0; i < matrix.size(); i++) {
            for (int j = 0; j < other.size().second; j++) {
                result.matrix[i][j] = 0; // Initialize the result element to 0
                for (int k = 0; k < matrix[0].size(); k++) {
                    result.matrix[i][j] += matrix[i][k] * other.matrix[k][j];
                }
            }
        }

        // Assign the result back to the current matrix
        matrix = result.matrix;

        return *this;
    }

    Matrix &operator*=(const T &scalar) {
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                matrix[i][j] *= scalar;
            }
        }
        return *this;
    }

    Matrix operator/(const T &scalar) const {
        if (scalar == 0.0) {
            throw std::invalid_argument("Cannot divide by zero");
        }
        Matrix result(matrix.size(), matrix[0].size());
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                result.matrix[i][j] = matrix[i][j] / scalar;
            }
        }
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
        return *this;
    }

    Matrix operator+(const Matrix &other) const {
        // Check if dimensions match
        if (matrix.size() != other.size().first || matrix[0].size() != other.size().second) {
            throw std::invalid_argument("Matrix dimensions don't match");
        }

        Matrix result(matrix.size(), matrix[0].size()); // Result matrix

        // Perform element-wise addition
        for (int i = 0; i < matrix.size(); i++) {
            for (int j = 0; j < matrix[0].size(); j++) {
                result.matrix[i][j] = matrix[i][j] + other.matrix[i][j];
            }
        }

        return result;
    }

    Matrix operator+(const T &scalar) const {
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
        for (int i = 0; i < matrix.size(); i++) {
            for (int j = 0; j < matrix[0].size(); j++) {
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

    Matrix operator-(const Matrix &other) const {
        // Check if dimensions match
        if (matrix.size() != other.size().first || matrix[0].size() != other.size().second) {
            throw std::invalid_argument("Matrix dimensions don't match");
        }

        Matrix result(matrix.size(), matrix[0].size()); // Result matrix

        // Perform element-wise subtraction
        for (int i = 0; i < matrix.size(); i++) {
            for (int j = 0; j < matrix[0].size(); j++) {
                result.matrix[i][j] = matrix[i][j] - other.matrix[i][j];
            }
        }

        return result;
    }

    Matrix operator-(const T &scalar) const {
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
        for (int i = 0; i < matrix.size(); i++) {
            for (int j = 0; j < matrix[0].size(); j++) {
                matrix[i][j] -= other.matrix[i][j];
            }
        }

        return *this;
    }

    Matrix &operator-=(const T &scalar) const {
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                matrix[i][j] -= scalar;
            }
        }
        return *this;
    }

    std::pair<size_t, size_t> size() const {
        return std::make_pair(matrix.size(), matrix[0].size());
    }

    T get(size_t row, size_t col) const {
        return matrix[row][col];
    }

    void set(size_t row, size_t col, T value) {
        matrix[row][col] = value;
    }
};
#endif
