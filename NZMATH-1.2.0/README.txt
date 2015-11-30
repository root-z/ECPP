NZMATH 1.2.0
=============

Introduction
------------

NZMATH is a Python based number theory oriented calculation system.
The centre of development in origin is Tokyo Metropolitan University.
Today it is developed at SourceForge.net.

This version 1.2.0 contains several new features, bug fixes.
The API can still be changed with versions.

Installation
------------

To install NZMATH on your computer, you must have Python 2.5 or
better.  If you don't have a copy of Python, please install it first.
Python is available from http://www.python.org/ .

The next step is to expand the NZMATH-1.2.0.tar.gz.  The way to do it
depends on your operating system.  On the systems with recent GNU tar,
you can do it with a single command::

 % tar xf NZMATH-1.2.0.tar.gz

where, % is the command line prompt.  Or with standard tar, you can do
it as::

 % gzip -cd NZMATH-1.2.0.tar.gz | tar xf -

Then, you have a child directory named NZMATH-1.2.0.

The third step is the last step, to install NZMATH to the standard
python path. Usually, this means to write files to somewhere under
/usr/lib or /usr/local/lib, and thus you have to have appropriate
write permission.  Typically, do as the following::

 % cd NZMATH-1.2.0
 % su
 # python setup.py install


Usage
-----

NZMATH is provided as a Python library package named 'nzmath', so
please use it as a usual package.  For more information please refer
Tutorial_.

.. _Tutorial: tutorial.html

Feedback
--------

Your feedbacks are always welcomed.  Please consider to join the
mailing list nzmath-user@tnt.math.se.tmu.ac.jp.  You can join the
list with writing a mail containing a line of "subscribe" in the body
to nzmath-user-request@tnt.math.se.tmu.ac.jp.  *Be careful* not to
send it to nzmath-user.


Copyright
---------
NZMATH is distributed under the BSD license.  See LICENSE.txt_ for
detail.

.. _LICENSE.txt: LICENSE.txt

Copyright (c) 2003-2012, NZMATH development group, all right reserved.
