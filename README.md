# NumericalMethods
A C++ console-based library for dealing with numerical methods.

## Contents
+ [About](#about)
+ [Installation](#installation)
+ [Usage](#usage)

## About <a name="about"></a>
The NumericalMethods is a library of `Matrix`, `Vector` and `SoLE` classes provided with features from the field of math of the same name. They are mostly designed as an easy-to-understand interfaces to be used in math calculations.

## Installation <a name="installation"></a>
There are several ways to install this library, depending on your development environment. Below provided instructions for users of Visual Studio and those who compile manually.

### Visual Studio
1. Copy the src folder to one of the directories on your PC.
2. Add .h and .cpp files as header and source files to your project.
3. Include the .h files in your source code using the `#include` directive, e.g.

```cpp
#include "Matrix.h"
#include "Vector.h"
#include "SoLE.h"
```

4. Compile and run your project.

### Manual
1. Copy the .h and .cpp files to the same directory as your main file.
2. Include the .h files in your source code, same as above for Visual Studio.
3. Compile the .cpp files using a C++ compiler (e.g., `g++`)

```bash
g++ -o your_executable your_main_file.cpp Matrix.cpp Vector.cpp SoLE.cpp
```

4. Run your executable.

Note: You should replace "your_executable" and "your_main_file" with your own executable and main file names.

## Usage <a name="usage"></a>

### Vector
```cpp
// main.cpp

// TASK: read two 4x1 vectors of integers from a files and output if they are equal (i.e., if they contain the same elements)

#include "Vector.h"


int main() {
    Vector<int> vec_1(4);
    vec_1.set_vector("vector_1.txt");

    Vector<int> vec_2(4);
    vec_2.set_vector("vector_2.txt");

    std::cout << vec_1 << std::endl;
    std::cout << vec_2 << std::endl;

    std::cout << vec_1 == vec_2 << std::endl;

    return 0;
}
```

### Matrix
```cpp
// main.cpp

// TASK: read a 3x3 matrix of integers from a file, output determinant and rank, transpose matrix and output the result

#include "Matrix.h"


int main() {
    Matrix<int> mat(3, 3);
    mat.set_matrix("matrix.txt");

    std::cout << mat.det() << std::endl;
    std::cout << mat.rank() << std::endl;

    mat.transpose();

    std::cout << mat << std::endl;

    // or in a single line: std::cout << mat.transposed() << std::endl;

    return 0;
}
```

### SoLE
```cpp
// main.cpp

// TASK: construct a <double> 5x5 SoLE with a matrix and a vector from files, solve the system (with any iterative method) and output the solution

#include "SoLE.h"


int main() {
    SoLE<double> sole(5, 5);
    sole.set_matrix("matrix.txt");
    sole.set_vector("vector.txt");

    Vector<int> solution = sole.solve_iter("seidel");

    std::cout << solution << std::endl;

    return 0;
}
```
