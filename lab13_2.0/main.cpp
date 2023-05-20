#include <iostream>
#include <fstream>
#include <vector>

class Matrix {
private:
    int rows;
    int cols;
    std::vector<std::vector<double>> data;

public:
    // Конструкторы
    Matrix() : rows(0), cols(0) {}
    Matrix(int rows, int cols) : rows(rows), cols(cols), data(rows, std::vector<double>(cols, 0.0)) {}
    Matrix(int rows, int cols, const std::vector<double>& data) : rows(rows), cols(cols), data(rows, std::vector<double>(cols, 0.0)) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                this->data[i][j] = data[i * cols + j];
            }
        }
    }



    // Методы создания матриц
    static Matrix zeros(int rows, int cols) {
        return Matrix(rows, cols); ч
    }

    static Matrix ones(int size) {
        Matrix mat(size, size);
        for (int i = 0; i < size; ++i) {
            mat(i, i) = 1.0;
        }
        return mat;
    }

    // Методы доступа к элементам матрицы
    double& operator()(int row, int col) {
        return data[row][col];
    }

    const double& operator()(int row, int col) const {
        return data[row][col];
    }

    // Операторы сравнения
    bool operator==(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            return false;
        }
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (data[i][j] != other.data[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    bool operator!=(const Matrix& other) const {
        return !(*this == other);
    }

    // Операторы математических операций
    Matrix operator+(const Matrix& other) const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] + other(i, j);
            }
        }
        return result;
    }

    Matrix operator-(const Matrix& other) const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] - other(i, j);
            }
        }
        return result;
    }

    Matrix operator*(double scalar) const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] * scalar;
            }
        }
        return result;
    }

    friend Matrix operator*(double scalar, const Matrix& mat) {
        return mat * scalar;
    }

    Matrix operator*(const Matrix& other) const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                double sum = 0.0;
                for (int k = 0; k < cols; ++k) {
                    sum += data[i][k] * other(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            rows = other.rows;
            cols = other.cols;
            data = other.data;
        }
        return *this;
    }

    // Методы ввода-вывода матриц
    friend std::istream& operator>>(std::istream& in, Matrix& mat) {
        int rows, cols;
        in >> rows >> cols;
        mat = Matrix(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                in >> mat(i, j);
            }
        }
        return in;
    }

    friend std::ostream& operator<<(std::ostream& out, const Matrix& mat) {
        out << mat.rows << "<-rows cols->" << mat.cols << std::endl << std::endl;
        for (int i = 0; i < mat.rows; ++i) {
            for (int j = 0; j < mat.cols; ++j) {
                out << mat(i, j) << " ";
            }
            out << std::endl;
        }
        return out;
    }

    // Методы чтения-записи матриц в файл
    void read(const std::string& filename) {
        std::ifstream infile(filename);
        infile >> *this;
    }

    void write(const std::string& filename) const {
        std::ofstream outfile(filename);
        outfile << *this;
    }
    //проверка на умножение
    bool is_multiplicable_with(const Matrix& other) const {
        return cols == other.rows;
    }
    //проверка на равенство размеров матрицы
    bool is_same_size_as(const Matrix& other) const {
        return (rows == other.rows) && (cols == other.cols);
    }
};

int main() {
    Matrix mat1(2, 3, { 1, 2, 3, 4,5,6 });
    Matrix mat2 = Matrix::ones(2);
    Matrix mat3;

    std::cin >> mat3;

    std::cout << "mat1: " << std::endl << mat1 << std::endl;
    std::cout << "mat2: " << std::endl << mat2 << std::endl;
    std::cout << "mat3: " << std::endl << mat3 << std::endl;



    if (mat1.is_multiplicable_with(mat2)) {
        Matrix mat4 = mat1 * mat2;
        std::cout << "mat4 = mat1 * mat2: " << std::endl << mat4 << std::endl;
    }
    else {
        std::cerr << "Matrices are not compatible for multiplication" << std::endl << std::endl;
    }


    if (mat1.is_same_size_as(mat3)) {
        Matrix mat5 = mat1 + mat3;
        std::cout << "mat6 = mat1 + mat3: " << std::endl << mat5 << std::endl;
    }
    else {
        std::cerr << "Matrices are not of the same size" << std::endl << std::endl;
    }



    if (mat2.is_same_size_as(mat3)) {
        Matrix mat6 = mat2 - mat3;
        std::cout << "mat6 = mat2 - mat3: " << std::endl << mat6 << std::endl;
    }
    else {
        std::cerr << "Matrices are not of the same size" << std::endl << std::endl;
    }


    Matrix mat7 = mat1 * 2;
    std::cout << "mat7 = mat1 * 2: " << std::endl << mat7 << std::endl;


    mat3.write("matrix.txt");
    Matrix mat8;
    mat8.read("matrix.txt");
    std::cout << "mat8 (from file): " << std::endl << mat8 << std::endl;

    return 0;
}