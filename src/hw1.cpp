#include "hw1.h"

namespace algebra 
{
    //zeros :
    Matrix zeros(size_t n, size_t m) {
        Matrix mat;
        for(size_t i{}; i<n; i++)
        {
            mat.push_back(std::vector<double>());
            for(size_t j{}; j<m; j++)
            {    
                mat[i].push_back(0);    
            }
        }
        return mat;
    }
    
    //ones :
    Matrix ones(size_t n, size_t m) {
        Matrix mat;
        for(size_t i{}; i<n; i++)
        {
            mat.push_back(std::vector<double>());
            for(size_t j{}; j<m; j++)
            {    
                mat[i].push_back(1);    
            }
        }
        return mat;
    }
   
    //random generator:
    Matrix random(size_t n, size_t m, double min, double max) {
        if(min > max){
            std::cout << "logic error, min cannot be bigger than max";
            throw std::logic_error("Error");
        }
        Matrix mat;
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(min, max);
        for(size_t i{}; i < n; i++){
            mat.push_back(std::vector<double>());
            for(size_t j{}; j < m; j++)
                mat[i].push_back(dist(mt));
        }
        return mat;
    }

    //show:
    void show(const Matrix& matrix){

        long unsigned int rows{ matrix.size() };            //number of rows in matrix
        long unsigned int columns{ matrix[0].size() };      //number of columns in matrix

        for (size_t i{}; i < rows; i++){
            for (size_t j{}; j < columns; j++){
                std::cout <<std::setw(3) << matrix[i][j];
            }
            std::cout << std::endl;
        }

    }

    //multiply
    Matrix multiply(const Matrix& matrix, double c){
        long unsigned int rows{ matrix.size() };            //number of rows in matrix
        long unsigned int columns{ matrix[0].size() };      //number of columns in matrix

        Matrix mat {matrix};
        for (size_t i{}; i < rows; i++)
            for(size_t j{}; j < columns; j++)
                mat[i][j] *= c;
            
        return mat;                
    }

    //multiply two matrices
    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2){
        Matrix mat1 { matrix1 };
        Matrix mat2 { matrix2 };
        
        if(mat1.empty() == true || mat2.empty() == true){
           /* if(mat1.empty())  
                return mat1;
            if(mat2.empty())
                return mat2;
            */
           Matrix mult{}; 
           return mult;        
        }      

        long unsigned int rows1{ mat1.size() };            //number of rows in matrix
        long unsigned int rows2{ mat2.size() };            //number of rows in matrix
        long unsigned int columns1{ mat1[0].size() };      //number of columns in matrix1
        long unsigned int columns2{ mat2[0].size() };      //number of columns in matrix2
        if (rows2 != columns1){
            std::cout << "You cannot multiply these two matrices";
            throw std::logic_error("Error");
            Matrix mult{};
            return mult;    
        }
        
        Matrix mult;                                       //multiply matrix


        for(int i{}; i < rows1; i++){
            mult.push_back(std::vector<double>());
            for(size_t j{}; j < columns2; j++)
                {
                    double temp {};
                    for(size_t k{}; k < columns1; k++){
                        temp += mat1[i][k] * mat2[k][j]; 
                    }
                    mult[i].push_back(temp);
                }        
        }        
        return mult;
    }

    //sum:
    Matrix sum(const Matrix& matrix, double c) {
        if(matrix.empty())
            return matrix;

        long unsigned int rows{ matrix.size() };            //number of rows in matrix
        long unsigned int columns{ matrix[0].size() };      //number of columns in matrix

        Matrix mat {matrix};
        for (size_t i{}; i < rows; i++)
            for(size_t j{}; j < columns; j++)
                mat[i][j] += c;
            
        return mat;                
    
    }

    //sum of two matrices:
    Matrix sum(const Matrix& matrix1, const Matrix& matrix2) {
        Matrix mat1 { matrix1 };
        Matrix mat2 { matrix2 };
        Matrix sum;                                       //multiply matrix
        if(mat1.empty() == true & mat2.empty() == true){
           Matrix mult{}; 
           return mult;        
        }
        if(mat1.empty() == true | mat2.empty() == true){
           throw std::logic_error("error");     
           Matrix mult{}; 
           return mult;        
        }      

        long unsigned int rows1{ mat1.size() };            //number of rows in matrix
        long unsigned int rows2{ mat2.size() };            //number of rows in matrix
        long unsigned int columns1{ mat1[0].size() };      //number of columns in matrix1
        long unsigned int columns2{ mat2[0].size() };      //number of columns in matrix2

        if (rows1 != rows2 || columns1 != columns2){
            std::cout << "You cannot sum these two matrices" << std::endl;
            throw std::logic_error("ERROR");
            return sum;    
        }
        sum = zeros(rows1, columns1);
        for(size_t i{}; i < rows1; i++)
            for(size_t j{}; j < columns1; j++)
                sum[i][j] = mat1[i][j] + mat2[i][j];
        
        return sum;
    }

    //transpose:
    Matrix transpose(const Matrix& matrix) {
        //Matrix mat {matrix};
        if(matrix.empty()){
            //throw std::logic_error("Error");
            return matrix;
        }

        Matrix trans;

        for(size_t i{}; i < matrix[0].size() ; i++) {
            trans.push_back(std::vector<double>());
            for(size_t j{}; j  < matrix.size(); j++ ){
                trans[i].push_back(matrix[j][i]);
            }
        }
        return trans;
    }

    //minor of matrix
    Matrix minor(const Matrix& matrix, size_t n, size_t m) {
        Matrix minor_matrix;
        unsigned int minor_row{}, minor_col{};
        if (matrix.empty() == true) {
            return matrix;
        }

        for(size_t k{}; k < matrix.size() - 1; k++){
            minor_matrix.push_back(std::vector<double>());
            for(size_t p{}; p < matrix[0].size() - 1; p++)
                minor_matrix[k].push_back(0);
        }

        for(size_t i{}; i < matrix.size(); i++){
            minor_row = i;
            if (i > n)
                minor_row--;
            
            for(size_t j{}; j < matrix[0].size(); j++) {
                minor_col = j;
                if(j > m)
                    minor_col--;
                
                if(i != n && j != m)
                    minor_matrix[minor_row][minor_col] = matrix[i][j];
            }
        }
        return minor_matrix;
    }

    //determinant of a matrix
    double determinant(const Matrix& matrix1){              //ye bug ee bood inja bekhater inke hame matrix arg o avaz nknm 
        Matrix matrix { matrix1 };                          //arg ro gozashtam matrix1 xD
        if(matrix.empty() == 1){
            double det { 1 };
            return det;
        }

        if(matrix.size() != matrix[0].size()){
            std::cout << "Matrix is not quadratic";
            throw std::logic_error("error");
        }
            
        double det {}; 
        if (matrix.size() == 1)
        {
            return matrix[0][0]; 
        }
        else if (matrix.size() == 2)
        {
            det = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
            return det;
        }
        else
        {
            for (int p = 0; p < matrix[0].size(); p++)
            {
                Matrix TempMatrix; 
                for (int i = 1; i < matrix.size(); i++)
                {
                    std::vector<double> TempRow;
                    for (int j = 0; j < matrix[i].size(); j++)
                    {
                        if (j != p)
                        {
                           TempRow.push_back(matrix[i][j]); 
                        }
                    }
                    if (TempRow.size() > 0)
                        TempMatrix.push_back(TempRow);
                }
                det = det + matrix[0][p] * pow(-1, p) * determinant(TempMatrix);
            }
            return det;
        }
    }
    
    //inverse:
    Matrix inverse(const Matrix& matrix){
        Matrix inv;
        Matrix invout;

        //empty matrix
        if (matrix.empty() == true) {
            //inv = {};
            return inv;
        }
        //0 determinant
        if (determinant(matrix) == 0){
            throw std::logic_error("error");
            inv={};
            return inv;

        }
        size_t rows = matrix. size();
        size_t cols = matrix[0]. size();
        
        //non squared matrices
        if (rows != cols){
            std::cout<<"this matrix does have an inverse"<<std::endl;
            throw std::logic_error ("Error");
           // inv = {};
            return matrix;
        }

        inv = zeros(rows,cols);
        invout=zeros(rows,cols);

        double sign { 1 };

        for (size_t i {}; i < rows; i++){
            for (size_t j {}; j < cols; j++){
                sign = (( i + j ) % 2 == 0)? 1: -1;
                inv[j][i] = sign * determinant(minor(matrix,i,j));      //darja transpose mikonim
                //std::cout << sign << std::endl;
            }
        //inv = transpose(inv);     
        }
        for (int i{}; i < rows; i++){
            for (int j{}; j < cols; j++){
            
                invout[i][j]=(inv[i][j] / determinant(matrix));
            }
        }
        
        return invout;    
        
    }    

    //concarenate
    Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis){
        Matrix conc_mat;
    /*  for(size_t i{}; i < mat1.size() + mat2.size(); i++){
            sum.push_back(std::vector<double>());
            for(size_t j{}; j < mat1[0].size(); j++){
                if (i < mat1.size()){
                    sum[i].push_back(mat1[i][j]);
                }
                else if(i >= mat1.size()){
                    sum[i].push_back(mat2[i][j]);
                }
            }                        
        } */
        if(axis == 0){
            if(matrix1[0].size() != matrix2[0].size()){
                std:: cout << "these two matrices don't have equal number of columns!";
                throw std::logic_error("Error");
                return conc_mat;
            }
            for(size_t i{}; i < matrix1.size() + matrix2.size(); i++){
                conc_mat.push_back(std::vector<double>());   
                for(size_t j{}; j < matrix1[0].size(); j++){
                    conc_mat[i].push_back(0);
                 }
            }
            unsigned long int n {matrix1.size()};
            for(size_t i{}; i < conc_mat.size(); i++){
                for(size_t j{}; j < conc_mat[0].size(); j++){
                    if(i < n)
                        conc_mat[i][j] = matrix1[i][j];
                    else if(i >= n)
                        conc_mat[i][j] = matrix2[i - n][j];
                }
            }
       }
        if(axis == 1){
            if(matrix1.size() != matrix2.size()){
                std::cout << "rows of two matrices are not equal!";
                throw std::logic_error("Error");
                return conc_mat;
            }

            for(size_t i{}; i < matrix1.size(); i++){
                conc_mat.push_back(std::vector<double>());   
                for(size_t j{}; j < matrix1[0].size() + matrix2[0].size(); j++){
                    conc_mat[i].push_back(0);
                }
            }
            unsigned long int n {matrix1[0].size()};
            for(size_t i{}; i < conc_mat.size(); i++){
                for(size_t j{}; j < conc_mat[0].size(); j++){
                    if(j < n)
                        conc_mat[i][j] = matrix1[i][j];
                    else if(j >= n)
                        conc_mat[i][j] = matrix2[i][j - n];
                }
            }
        }
        return conc_mat;
    }
    
    //swap rows
    Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2){
        if(matrix.empty()){
            return matrix;
        }
        
        if (matrix.size() < r1 + 1 || matrix.size() < r2 + 1)
        {
            std::cout << " check your matrix dimension with input rows " << std::endl;
            throw std::logic_error("Error");
        }
        
        Matrix swaped {matrix};
        for(size_t i{}; i < matrix[0].size(); i++){
            double temp{ swaped[r1][i] };
            swaped[r1][i] = swaped[r2][i];
            swaped[r2][i] = temp;
        }
        return swaped;
    }

    //multiply a row to a number
    Matrix ero_multiply(const Matrix& matrix, size_t r, double c){
        if(matrix.empty()){
            return matrix;
        }
        Matrix mult { matrix };
        for(size_t i{}; i < mult[0].size(); i++){
            mult[r][i] *= c;
        }
        return mult;
    }
    
    //combination of sum an multiply
    Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2){
        if(matrix.empty()){
            return matrix;
        }

        if (matrix.size() < r1 + 1 || matrix.size() < r2 + 1){
            std::cout << " check your matrix dimension with input rows " << std::endl;
            throw std::logic_error("Error");
        }

        Matrix sum {matrix};        //khode matrix ro rikhtam tosh ta dimantion ash yeksan bashan
        for(size_t i{}; i < sum[0].size(); i++)
            sum[r2][i] = matrix[r1][i] * c + matrix[r2][i];
        return sum;
    }

    //make upper triangle
    Matrix upper_triangular(const Matrix& matrix){
        Matrix trimat { matrix };
        if(matrix.empty()){
            return matrix;
        }
        if(matrix.size() != matrix[0].size()){
            std::cout << "non square matrices cannot become upper triangle" << std::endl;
            throw std::logic_error("Error");
        }

        for(size_t i{}; i < trimat[0].size(); i++){
            for(size_t j{ i + 1 }; j < trimat.size(); j++){
                if(i + j > trimat.size())
                    break;
                double temp {-1 * trimat[j][i] / trimat[i][i]};
                trimat = ero_sum(trimat, i, temp, j);
            }
        }
        return trimat;
    }
}
