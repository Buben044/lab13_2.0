#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

template<typename T>
class Matrix {
private:
    int rows;
    int cols;
    std::vector<std::vector<T>> data;

public:
    Matrix<T>() : rows(0), cols(0) {}
    Matrix<T>(int rows, int cols) : rows(rows), cols(cols), data(rows, std::vector<T>(cols, 0.0)) {}
    Matrix<T>(int rows, int cols, const std::vector<T>& data) : rows(rows), cols(cols), data(rows, std::vector<T>(cols, 0.0)) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                this->data[i][j] = data[i * cols + j];
            }
        }
    }
    
    void insert(T n, int i, int j) {
        if (rows < i + 1) {
            data.resize(rows + 1);
            rows++;
        }if (data[i].size() < j + 1) data[i].push_back(n);
        else data[i][j] = n;
        if (cols < j + 1) cols++;
    }
    
    static Matrix<T> zeros(int rows, int cols) {
        return Matrix<T>(rows, cols);
    }
    
    static Matrix<T> ones(int size) {
        Matrix<T> mat(size, size);
        for (int i = 0; i < size; ++i) {
            mat(i, i) = (T) 1;
        }
        return mat;
    
    
    }T& operator()(int row, int col) {
        return data[row][col];
    }
    
    const T& operator()(int row, int col) const {
        return data[row][col];
    }
    
    bool operator==(const Matrix<T>& other) const {
        if (rows != other.rows || cols != other.cols) {
            return false;
        }for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (data[i][j] != other.data[i][j]) {
                    return false;
                }
            }
        }return true;
    }
    bool operator!=(const Matrix<T>& other) const {
        return !(*this == other);
    }
    
    Matrix<T> operator+(const Matrix<T>& other) const {
        if (!(this->is_same_size_as(other))) {
            std::cerr << "Matrices are not of the same size" << std::endl;
            return *this;
        }
        Matrix<T> result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] +  (T) other(i, j);
            }
        }return result;
    }
    
    Matrix<T> operator-(const Matrix<T>& other) const {
        if (!(this->is_same_size_as(other))) {
            std::cerr << "Matrices are not of the same size" << std::endl;
            return *this;
        }Matrix<T> result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] - (T) other(i, j);
            }
        }return result;
    }
    
    Matrix<T> operator*(double scalar) const {
        Matrix<T> result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = (double)data[i][j] * scalar;
            }
        }return result;
    }
    
    friend Matrix<T> operator*(double scalar, const Matrix<T>& mat) {
        return mat * scalar;

    }
    
    Matrix<T> operator*(const Matrix<T>& other) const {
        if (!(this->is_multiplicable_with(other))) {
            std::cerr << "Matrices are not compatible for multiplication" << std::endl;
            return *this;
        }
        Matrix<T> result(rows, other.cols);
        for (int i = 0; i < rows; ++i) {
            for (int k = 0; k < other.cols; ++k) {
                T sum = 0;
                for (int j = 0; j < cols; ++j) {
                    sum += data[i][j] * other(j, k);
                }
                result(i, k) = sum;
            }
        }return result;
    }
    
    Matrix<T>& operator=(const Matrix<T>& other) {
        if (this != &other) {
            rows = other.rows;
            cols = other.cols;
            data = other.data;
        }return *this;
    }
    
    friend std::istream& operator>>(std::istream& in, Matrix<T>& mat) {
        int rows, cols;
        in >> rows >> cols;
        mat = Matrix<T>(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                in >> mat(i, j);
            }
        }return in;
    }
    
    friend std::ostream& operator<<(std::ostream& out, const Matrix<T>& mat) {
        out << mat.rows << " " << mat.cols << std::endl << std::endl;
        for (int i = 0; i < mat.rows; ++i) {
            for (int j = 0; j < mat.cols; ++j) {
                out << mat(i, j) << " ";
            }
            out << std::endl;
        }
        return out;
    }
    
    void read(const std::string& filename) {
        std::ifstream infile(filename);
        infile >> *this;
    }
    
    void write(const std::string& filename) const {
        std::ofstream outfile(filename);
        outfile << *this;
    }
    
    bool is_multiplicable_with(const Matrix<T>& other) const {
        return cols == other.rows;
    }
    
    bool is_same_size_as(const Matrix<T>& other) const {
        return (rows == other.rows) && (cols == other.cols);
    }
    
    double det() {
        if (rows != cols) {
            std::cout << "The Matrix<T> must be square" << std::endl << std::endl;
            return 0;
        }if (rows == 1) return (double)data[0][0];
        if (rows == 2) return data[0][0] * data[1][1] - data[0][1] * data[1][0];
        double res = 0;
        for (int i = 0; i < rows; i++) {
            double temp = 1;
            for (int j = 0; j < rows; j++) {
                temp *= (double)data[j][(i + j) % rows];
            }res += temp;
        }for (int i = 0; i < rows; i++) {
            double temp = 1;
            for (int j = 0; j < rows; j++) {
                temp *= (double)data[j][(rows + i - j) % rows];
            }res -= temp;
        }return res;
    }
    
    Matrix<T> get_small(int i, int j) {
        Matrix<T> m;
        bool flag = false;
        for (int k = 0; k < rows; k++) {
            if (k != i) {
                bool flag1 = false;
                for (int g = 0; g < cols; g++) {
                    if (g != j) {
                        int k1, g1;
                        if (flag) k1 = k - 1;
                        else k1 = k;
                        if (flag1) g1 = g - 1;
                        else g1 = g;
                        m.insert(data[k][g], k1, g1);
                    }
                    else flag1 = true;
                }
            }
            else flag = true;
        }return m;
    }
    
    Matrix<T> transpon() {
        Matrix<T> m;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                m.insert(data[i][j], j, i);
            }
        }return m;
    }
    
    Matrix<T> findAlgebraicComplement() {
        if (rows != cols) {
            std::cout << "The Matrix<T> must be square" << std::endl << std::endl;
            return *this;
        }Matrix<T> m;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                m.insert(this->get_small(i, j).det() * pow((-1), i + j), i, j);
            }
        }return m;
    }
    
    Matrix<T> operator!() {
        if (!(this->det())) {
            std::cout << "The Matrix<T> must be square" << std::endl << std::endl;
            return *this;
        }if (this->is_E()) return *this;
        if (this->is_O()) return *this;
        Matrix<T> m = (this->findAlgebraicComplement().transpon() / this->det());
        return m;
    }
};

int main() {

    Matrix<double> mat1(2, 3, { 1, 2, 3, 4,5,6 });
    Matrix<double> mat2 = Matrix<double>::ones(2);
    Matrix<double> mat3;

    std::cin >> mat3;

    std::cout << "mat1: " << std::endl << mat1 << std::endl;
    std::cout << "mat2: " << std::endl << mat2 << std::endl;
    std::cout << "mat3: " << std::endl << mat3 << std::endl;


    Matrix<double> mat4 = mat1 * mat2;
    std::cout << "mat4 = mat1 * mat2: " << std::endl << mat4 << std::endl;

    Matrix<double> mat5 = mat1 + mat3;
    std::cout << "mat6 = mat1 + mat3: " << std::endl << mat5 << std::endl;

    Matrix<double> mat6 = mat2 - mat3;
    std::cout << "mat6 = mat2 - mat3: " << std::endl << mat6 << std::endl;

    
    Matrix<double> mat7 = mat1 * 2.0;
    std::cout << "mat7 = mat1 * 2: " << std::endl << mat7 << std::endl;


    mat3.write("Matrix.txt");
    Matrix<int> mat8;

    mat8.read("Matrix.txt");
    std::cout << "mat8 (from file): " << std::endl << mat8 << std::endl;

    Matrix<double> mat10 = mat1;
    mat10 = mat10.transpon();

    Matrix<double> mat9 = mat10 * mat1;
    std::cout << "mat9 = mat1 * mat10: " << std::endl << mat9 << std::endl;


    return 0;
}
