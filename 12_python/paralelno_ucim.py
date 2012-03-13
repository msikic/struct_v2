#!/usr/bin/env python

from multiprocessing import Pool

def f(x):
	return x*x

if __name__ == '__main__':
	pool = Pool(processes=4)                  
	print pool.map(f, range(50)), pool.pid

