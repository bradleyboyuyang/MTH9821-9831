## Code structure

- `q1_writer.cpp`: question 1, write a shared integer into the memory
- `q1_reader.cpp`: question 1, read the shared integer from the memory and release the memory
- other files: belong to question 2 and question 3



## How to run

1. Use `CMakeLists.txt` to generate makefile, build the project and three executable files will be generated (`writer`, `reader`, `test_products`) under folder `build`.
2. First run `writer`, then run `reader`, this shows the solution of question 1.
3. Run `test_products`, it contains the testing output for question 2 (retrieve three futures from the `futureProductService`) and question 3 (functionality of `getBond()` and `getFuture()`)