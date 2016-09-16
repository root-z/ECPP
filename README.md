# ECPP

To run the program, two libraries: nzmath and mpmath are required.
If you don't have the libraries installed, in the project folder run the following command:

1. pip install virtualenv

2. virtualenv -p /usr/bin/python2.7 venv

3. source venv/bin/activate

4. pip install mpmath

5. pip install nzmath

Note this process is for Linux/OSX. For windows virtualenv requires slight different command to activate and pip might not work for nzmath.
If Windows machine is used please go to the website of these libraries and install as suggested. The script has been tested on Windows.

After installation of the program, there are two ways to run the program.

1. In Bash/Terminal/..., Enter
    python prime_test.py (prime_number)

2. In Interactive Python or scripts
    from prime_test import prime
    prime(prime_number)

In the latest version I have also added a Flask web service for running the prime tester. See ecpp_service.py.