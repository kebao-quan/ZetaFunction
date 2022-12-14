#to use:
#first, get a list of primes.
#for example, primeList=p.primelist(1000000)
#
#second, call NPioverL on various large values of l
#for example, NPioverL(100000,primeList,2)

import math
import numpy as np
import time
import os
from decimal import *
from multiprocessing import Pool, RawArray
# def primelist(max):
# 	pl=[2]
# 	i=3
# 	while(i<max):
# 		j=2
# 		while(j<=math.sqrt(i)):
# 			if (i%j==0):
# 				break
# 			j=j+1
# 		if j>math.sqrt(i):
# 			pl=pl+[i]
# 		i=i+2
# 	return pl

def primelist(max):
    #Input max>=6, Returns a array of primes, 2 <= p < n
    sieve = np.ones(max//3 + (max%6==2), dtype=bool)
    for i in range(1,int(max**0.5)//3+1):
        if sieve[i]:
            k=3*i+1|1
            sieve[k*k//3::2*k] = False
            sieve[k*(k-2*(i&1)+4)//3::2*k] = False
    return np.r_[2,3,((3*np.nonzero(sieve)[0][1:]+1)|1)]


def NPioverL(l, power):
    return N(l,power)*(math.pi**power)/(l**power)

Prime_dict = {}
def init_worker(primeListRaw_, X_shape):
    Prime_dict["X"] = primeListRaw_
    Prime_dict["X_shape"] = X_shape
    
    

def N(l,power):
	maxterm = 15
	if l <= 10**10:
		maxterm = 12
	if l <= 10**9:
		maxterm = 11
	if l <= 10**8:
		maxterm = 9 #correct
	if l <= 10**7:
		maxterm = 9 #correct
	term=1
	sign=-1
	sum=l**power

	function_args = []
	while term < maxterm:
		function_args.append((l,term,power))
		term+=1
 
	with Pool(processes=os.cpu_count(),initializer=init_worker, initargs=(PrimeListRaw_a, X_shape)) as p:
		rets = p.map(Nterm, function_args)
		p.close()

	rets.sort(reverse=True)
	for value in rets:
		sum = sum + value * sign
		sign = sign * (-1)
	print("Terms: ", rets)
	return sum

def Nterm(args):
    
	X_np = np.frombuffer(Prime_dict['X']).reshape(Prime_dict['X_shape'])
	prime_list = X_np.ravel()
    
	l,term,power = args
	del args
	primeIndices=[-1]*term
	ret = NtermRecursive(l,primeIndices,prime_list,power,term)
	return ret

def NtermRecursive(l,primeIndices,primeList,power,term):
	i=0
	while i<len(primeIndices) and primeIndices[i]!=-1:
		i+=1

	if i==len(primeIndices):
		return (math.floor(l/product(primeIndices,primeList)))**power

	sum=0
	if i==0 and primelimit >= 1:
		while True:
			primeIndices[0]+=1

			if primeList[primeIndices[i]] > (l/leastPrimeProduct(term-1,primeList)):
				break

			sum+=NtermRecursive(l,[*primeIndices],primeList,power,term)
		return sum

	if i==1 and primelimit >= 2:
		while True:
			primeIndices[1]+=1

			if primeIndices[1]>=primeIndices[0]:
				break

			if primeList[primeIndices[1]]*primeList[primeIndices[0]] > (l/leastPrimeProduct(term-2,primeList)):
				break

			sum+=NtermRecursive(l,[*primeIndices],primeList,power,term)
		return sum

	if i==2 and primelimit >= 3:
		while True:
			primeIndices[2]+=1

			if primeIndices[2]>=primeIndices[1]:
				break

			if primeList[primeIndices[2]]*primeList[primeIndices[1]]*primeList[primeIndices[0]] > (l/leastPrimeProduct(term-3,primeList)):
				break

			sum+=NtermRecursive(l,[*primeIndices],primeList,power,term)
		return sum

	if i==3 and primelimit >= 4:
		while True:
			primeIndices[3]+=1

			if primeIndices[3]>=primeIndices[2]:
				break

			if primeList[primeIndices[3]]*primeList[primeIndices[2]]*primeList[primeIndices[1]]*primeList[primeIndices[0]] > (l/leastPrimeProduct(term-4,primeList)):
				break

			sum+=NtermRecursive(l,[*primeIndices],primeList,power,term)
		return sum

	if i==4 and primelimit >= 5:
		while True:
			primeIndices[4]+=1

			if primeIndices[4]>=primeIndices[3]:
				break

			if primeList[primeIndices[4]]*primeList[primeIndices[3]]*primeList[primeIndices[2]]*primeList[primeIndices[1]]*primeList[primeIndices[0]] > (l/leastPrimeProduct(term-5,primeList)):
				break

			sum+=NtermRecursive(l,[*primeIndices],primeList,power,term)
		return sum

	while True:
		primeIndices[i]+=1

		if primeIndices[i]>=primeIndices[i-1]:
			break

		if product(primeIndices[:i+1],primeList)>l:
			break

		sum+=NtermRecursive(l,[*primeIndices],primeList,power,term)
	return sum

def leastPrimeProduct(n,primeList):
    p=1
    for i in range(n):
        p=p*primeList[i]
    return p

def product(primeIndices,primeList):
	p=1
	for primeIndice in primeIndices:
		p=p*primeList[primeIndice]
	return p










L = 10**8
power = 2


primeList_ = primelist(L+10000)
X_shape = (1, len(primeList_))
primeList_1 = np.reshape(primeList_, X_shape)
PrimeListRaw_a = RawArray('Q', len(primeList_))

X_np = np.frombuffer(PrimeListRaw_a).reshape(X_shape)
np.copyto(X_np, primeList_1)

primelimit = 2					#the number of prime to restrict
# loadPrimeList = False 			#set this to True if you want to use existing primeList
# StorePrimeLocal = False			#set this to True if you want to save the primeList

if __name__ == '__main__':
# 	print("maximum threads: ", os.cpu_count())
# 	if loadPrimeList:
# 		with open(r'PrimeList.txt', 'r') as fp:
# 			primeList_ = list(map(int, fp.readline().split()))
# 		print("read PrimeList.txt complete")
# 	else:
# 		start_time = time.time()
# 		print("--- %s seconds --- for getting the prime list" % (time.time() - start_time))
# 		print("len of primeList: ", len(primeList_))
# 		print("------------------------------------------------\n")
# 		if StorePrimeLocal:
# 			with open(r'PrimeList.txt', 'w') as fp:
# 				for item in primeList_:
# 					fp.write("%s " % item)

	start_time = time.time()
	print("timer start! Computing L =", L, "  ---  power =", power, "  ---  Limit prime =", primelimit, "......")
	result = NPioverL(L, power)
	print("result: ", Decimal(result))
	print("--- %s seconds --- for getting the result" % (time.time() - start_time))
	print("timer end!")