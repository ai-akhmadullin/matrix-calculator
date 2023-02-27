import matrix_calc
from typing import List


def output_description(operation_count: int) -> None:
    """Print available matrix operations."""
    print("Matrix operations currently available:\n")
    for i in range(operation_count):
        print("{} - {}".format(i + 1, matrix_calc.operations[i]))


def get_operation(operation_count: int) -> int:
    """Prompt the user to enter the desired operation and return its number."""
    answer = 1
    while True:
        try:
            answer = int(input("\nEnter the number of the operation you want to do -> "))
            if answer > operation_count or answer < 1:
                raise ValueError

        # If the number entered by the user is out of range of possible operations.
        except ValueError:
            print("\nSuch option is not found, please try again")
            continue
        else:
            break
    print("You chose the following operation: {}".format(matrix_calc.operations[answer - 1]))
    return answer


def input_matrix(matrix_num: int) -> List[List[List[float]]]:
    """Prompt the user to enter matrices and return a list of entered matrices."""
    print("\nMatrix ends with an empty string")
    print("Start entering the matrix:")

    # List with all entered matrices.
    tensor = []
    for i in range(matrix_num):
        if i != 0:
            print("Next matrix:")
        end = False

        # The number of columns must be equal in a matrix and we will check for that. Since the matrix
        # has not been entered yet, we will set 'cols' to be -1 at the beginning and change it later.
        cols = -1
        _matrix = []
        while not end:
            row = input()

            # An empty string means the end of the matrix input.
            if row == '':
                end = True
            else:
                row = [float(x) for x in row.split()]

                # If the matrix has the wrong size.
                if (len(row) != cols) and (cols != -1):
                    raise ValueError("The number of elements in rows must be equal. Start again.")

                # In case of successful row entry, the size of the last row = the size of the current row.
                cols = len(row)
                _matrix.append(row)
        tensor.append(_matrix)

    return tensor


def apply_operation(operation: int) -> None:
    """Apply the operation and output the result."""
    operation_func = [matrix_calc.addition,
                      matrix_calc.subtraction,
                      matrix_calc.multiply,
                      matrix_calc.scalar,
                      matrix_calc.transpose,
                      matrix_calc.invert,
                      matrix_calc.determinant,
                      matrix_calc.lu,
                      matrix_calc.plu,
                      matrix_calc.rref]

    matrices = 1

    # The first three operations can be performed with several matrices.
    if operation <= 3:
        matrices = int(input("\nHow many matrices are you going to operate on? -> "))

    # Prompt the user to enter matrices.
    tensor = input_matrix(matrices)

    # Perform the operation.
    result = operation_func[operation - 1](tensor)

    print("Result:")

    # The determinant operation returns a number, so we just print it.
    if operation == 7:
        print(result)

    # The decomposition operation returns multiple matrices, so we print them one by one.
    elif operation == 8:
        for i, j in enumerate(['L', 'U']):
            print("{} = ".format(j))
            result[i].print()
    elif operation == 9:
        for i, j in enumerate(['P', 'L', 'U']):
            print("{} = ".format(j))
            result[i].print()

    # In all other cases, the output is a matrix, which is printed by the Matrix.print method.
    else:
        result.print()


def update_operation(message: str, operation: int, operation_count: int) -> int:
    print(message)
    answer = input("\nWould you like to change the operation? [y/n] ")
    if answer.lower() == 'y':
        operation = get_operation(operation_count)
    return operation


def main() -> None:
    """Carry out the main stages of the program."""
    # Display possible operations with matrices.
    operation_count = len(matrix_calc.operations)
    output_description(operation_count)

    # Prompt the user to enter the number of the operation.
    operation = get_operation(operation_count)

    # Apply selected operation.
    while True:
        try:
            apply_operation(operation)
        except ZeroDivisionError:
            message = "Division by zero, LU decomposition cannot be carried out. Use PLU decomposition instead or " \
                      "type in another matrix. "
            operation = update_operation(message, operation, operation_count)
            continue
        except ValueError as e:
            operation = update_operation(e, operation, operation_count)
            continue
        else:
            break


main()
