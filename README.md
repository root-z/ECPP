# ECPP

To run the program, two libraries: nzmath and mpmath are required.
If you don't have the libraries installed, in the project folder run the following command:
pip install virtualenv
virtualenv venv
source venv/bin/activate
pip install mpmath
pip install nzmath

Note this process is for Linux/OSX. For windows virtualenv requires slight different command to activate and pip might not work for nzmath.
If Windows machine is used please go to the website of these libraries and install as suggested. The script has been tested on Windows.

After installation of the program, there are two ways to run the program.

1. In Bash/Terminal/..., Enter
    python prime_test.py (prime_number)

2. In Interactive Python or scripts
    from prime_test import prime
    prime(prime_number)