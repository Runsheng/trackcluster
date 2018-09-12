# -*- coding: utf-8 -*-
"""
@package    ssw
@brief      **SSW** is a python interface with a higly efficient Smith-Waterman library written in C language. A dynamic library add to be compliled
to interface C and python by using the Makefile. This will generate libssw.so require for proper program execution.

* The python wrapper 'ssw_wrapp" allow to perform high performance pairwise DNA alignment and return a simple python object describing the alignement.
* pyssw.py is a standalone high level interface to send multiple queries to a reference than return results as a sam file
* C sources ssw.c and ssw.h were forked and modified from Mengyao's [Complete-Striped-Smith-Waterman-Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)

@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

__all__ = ["ssw_wrap", "pyssw"]
