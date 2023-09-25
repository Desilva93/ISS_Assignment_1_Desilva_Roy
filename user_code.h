#ifndef USER_CODE_H
#define USER_CODE_H

#include <vector>
#include <string>
#include <utility>
#include <unordered_map>
#include <stdio.h>

// Function to convert a dense matrix to a sparse matrix representation
void denseToSparse(const std::vector<std::vector<int> >& denseMatrix,
                   std::vector<int>& values,
                   std::vector<int>& columnIndices,
                   std::vector<int>& rowPointers){
                    // Write your code here
                    // Identify all of the non-zero elements in the dense matrix.
                    std::vector<std::pair<int, int>> non_zero_elements;
                    for (int i = 0; i < denseMatrix.size(); i++) {
                        for (int j = 0; j < denseMatrix[i].size(); j++) {
                        if (denseMatrix[i][j] != 0) {
                            non_zero_elements.push_back({i, j});
                        }
                        }
                    }

                    // Create the CSR matrix.
                    values.resize(non_zero_elements.size());
                    columnIndices.resize(non_zero_elements.size());
                    rowPointers.resize(denseMatrix.size() + 1);

                    // Count the number of non-zero elements in each row.
                    for (int i = 0; i < non_zero_elements.size(); i++) {
                        rowPointers[non_zero_elements[i].first + 1]++;
                        }

                    // Calculate the row pointers.
                    for (int i = 1; i < rowPointers.size(); i++) {
                        rowPointers[i] += rowPointers[i - 1];
                        }      

                    // Fill in the CSR matrix.
                    for (int i = 0; i < non_zero_elements.size(); i++) {
                        values[rowPointers[non_zero_elements[i].first] + i - rowPointers[non_zero_elements[i].first]] = denseMatrix[non_zero_elements[i].first][non_zero_elements[i].second];
                        columnIndices[rowPointers[non_zero_elements[i].first] + i - rowPointers[non_zero_elements[i].first]] = non_zero_elements[i].second;
                        }
                    
                    return;
}


//Function to perform dense matrix multiplication

void denseMatrixMultiplication(const std::vector<std::vector<int> >& denseMatrix1,
								const std::vector<std::vector<int> >& denseMatrix2,
                                int squareMatrixSize,
                                std::vector<std::vector<int> >& denseMatrix_result){

            for (int i = 0; i < squareMatrixSize; i++) {
                for (int j = 0; j < squareMatrixSize; j++) {
                    denseMatrix_result[i][j] = 0;
                    for (int k = 0; k < squareMatrixSize; k++) {
                        denseMatrix_result[i][j] += denseMatrix1[i][k] * denseMatrix2[k][j];
                    }
                }
            }
			return;
}


//Function to perform dense matrix multiplication using block matrix multiplication
void blockMatrixMultiplication(const std::vector<std::vector<int> >& denseMatrix1,
								const std::vector<std::vector<int> >& denseMatrix2,
                                int N, int B,
                                std::vector<std::vector<int> >& denseMatrix_result){
                                // Write your code here
                                for (int i = 0; i < N; i += B) {
                                    for (int j = 0; j < N; j += B) {
                                    // Iterate over the blocks of the second matrix.
                                        for (int k = 0; k < N; k += B) {
                                            // Calculate the block element at row i, column j of the result matrix.
                                            for (int ii = i; ii < i + B; ii++) {
                                                for (int jj = j; jj < j + B; jj++) {
                                                    for (int kk = k; kk < k + B; kk++) {
                                                        denseMatrix_result[ii][jj] += denseMatrix1[ii][kk] * denseMatrix2[kk][jj];
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }                                
                                return;
                            }

// Function to perform sparse matrix multiplication
void sparseMatrixMultiplication(const std::vector<int>& values1,
                                const std::vector<int>& columnIndices1,
                                const std::vector<int>& rowPointers1,
                                const std::vector<int>& values2,
                                const std::vector<int>& columnIndices2,
                                const std::vector<int>& rowPointers2,
                                int N,
                                std::vector<int>& valuesResult,
                                std::vector<int>& columnIndicesResult,
                                std::vector<int>& rowPointersResult){
                                // Write your code here
                                // Initialize the result matrix.
                                valuesResult.clear();
                                columnIndicesResult.clear();
                                rowPointersResult.resize(N + 1, 0);

                                // Iterate over the rows of the first matrix.
                                for (int i = 0; i < N; i++) {
                                    // Iterate over the elements in the current row of the first matrix.
                                    for (int j = rowPointers1[i]; j < rowPointers1[i + 1]; j++) {
                                    // Get the column index of the current element.
                                        int column_index1 = columnIndices1[j];

                                        // Iterate over the elements in the corresponding column of the second matrix.
                                        for (int k = rowPointers2[column_index1]; k < rowPointers2[column_index1 + 1]; k++) {
                                            // Get the column index of the current element.
                                            int column_index2 = columnIndices2[k];

                                            // Calculate the product of the two elements.
                                            int product = values1[j] * values2[k];

                                            // If the product is non-zero, add it to the result matrix.
                                            if (product != 0) {
                                            valuesResult.push_back(product);
                                            columnIndicesResult.push_back(column_index2);
                                            }
                                        }
                                    }
                                }

                                // Calculate the row pointers for the result matrix.
                                for (int i = 1; i < rowPointersResult.size(); i++) {
                                    rowPointersResult[i] = rowPointersResult[i - 1] + valuesResult.size();
                                }                                
                                return;
}


//function For finding patterns in graphs using a map
void findTrianglesAsMap(const std::vector<std::vector<int> >& adjMat, int N, std::map<std::vector<int>, int>& outputMap){
    // Write your code here
    // Iterate over all the rows of the adjacency matrix.
    for (int i = 0; i < N; i++) {
    // Iterate over all the columns of the adjacency matrix.
        for (int j = i + 1; j < N; j++) {
            // If there is an edge between the current row and column, then iterate over all the columns again.
            if (adjMat[i][j] != 0) {
            for (int k = j + 1; k < N; k++) {
                // If there is an edge between the current row and column, and there is an edge between the current column and the previous column, then we have found a triangle.
                if (adjMat[i][k] != 0 && adjMat[j][k] != 0) {
                // Sort the edges of the triangle in ascending order.
                std::vector<int> edges = {adjMat[i][j], adjMat[j][k], adjMat[i][k]};
                std::sort(edges.begin(), edges.end());

                // Increment the count of the triangle pattern.
                outputMap[edges]++;
                }
            }
            }
        }
    }    
    return;
}


//function For finding patterns in graphs using a vector

void findTrianglesAsVec(const std::vector<std::vector<int> >& adjMat, int N,  std::vector<std::vector<int> >& outputVec){
    // Write your code here
    // Initialize a hash table to store the triangle patterns and their counts.
    std::map<std::vector<int>, int> triangleCounts;

    // Iterate over all the rows of the adjacency matrix.
    for (int i = 0; i < N; i++) {
        // Iterate over all the columns of the adjacency matrix.
        for (int j = i + 1; j < N; j++) {
        // If there is an edge between the currmkdirent row and column, then iterate over all the columns again.
        if (adjMat[i][j] != 0) {
            for (int k = j + 1; k < N; k++) {
            // If there is an edge between the current row and column, and there is an edge between the current column and the previous column, then we have found a triangle.
            if (adjMat[i][k] != 0 && adjMat[j][k] != 0) {
                // Sort the edges of the triangle in ascending order.
                std::vector<int> edges = {adjMat[i][j], adjMat[j][k], adjMat[i][k]};
                std::sort(edges.begin(), edges.end());

                // Increment the count of the triangle pattern.
                triangleCounts[edges]++;
            }
            }
        }
        }
    }

    // Iterate over the triangle patterns and their counts, and add them to the output vector.
    for (const auto& [pattern, count] : triangleCounts) {
        outputVec.push_back({pattern[0], pattern[1], pattern[2], count});
    }
    return;
}

#endif // USER_CODE_H
