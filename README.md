# Python Project Template
A python project template that uses tox to run test suites. setup.cfg lists metadata for the package and lists the required packages for the program. It also lists the requirements for the testing environment. 
Tox.ini creates all of the environments (right now it's just for py38 and testenv). testenv inside of tox.ini tells tox to look in setup.cfg to find dependencies for testing and also gives the command line instruction for testing.

For more information, see [this link](https://packaging-guide.openastronomy.org/en/latest/index.html)
