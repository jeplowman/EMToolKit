# examples/add_two.py

"""
Adding Two Numbers
==================

This example demonstrates how to use the add_numbers function from EMToolKit.

The function takes two numbers as input and returns their sum.
"""

import matplotlib.pyplot as plt
from numpy import add
# from EMToolKit.example_module import add_numbers

def main():
    # Example usage of the add_numbers function

    print("Hello World")

    result1 = add(2, 3)
    result2 = add(4.5, 5.5)

    # Print the results
    print(f"add_numbers(2, 3) = {result1}")
    print(f"add_numbers(4.5, 5.5) = {result2}")

    # Plot the results
    plt.figure()
    plt.bar(["2 + 3", "4.5 + 5.5"], [result1, result2])
    plt.title("Results of add_numbers function")
    plt.xlabel("Operations")
    plt.ylabel("Results")
    plt.show()

if __name__ == "__main__":
    main()