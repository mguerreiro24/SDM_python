#############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++Built by Miguel Fernandes Guerreiro+++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++22/05/2017++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#############################################################################
'''
sampling iterable objects
'''
from random import randint
from copy import deepcopy
def WithReposition(n,iterable):
   '''subsamples **with** reposition
requires:[iterable: str(), tuple(), list()] and [n: int()] number of samples
to retrieve
ensures:length n list with n random samples from iterable provided
'''
   l = len(iterable)
   k = []
   while n>len(k):
      k.append(randint(0, l-1))

   l = []
   wk = deepcopy(iterable)
   for i in list(k):
      l.append(wk[i])
   return l


def WithoutReposition(n,iterable):
   '''subsamples **without** reposition
requires:[iterable: str(), tuple(), list()] and [n: int()] number of samples
to retrieve
ensures:length n list with n different random samples from iterable provided
'''
   l = len(iterable)
   k = set()
   while n>len(k):
      k.add(randint(0, l-1))
   l = []
   wk = deepcopy(iterable)
   for i in list(k):
      l.append(wk[i])
   return l


def WithoutReposition2(n,iterable):
   '''subsamples **without** reposition
requires:[iterable: str(), tuple(), list()] and [n: int()] number of samples
to retrieve
ensures:length n list with n different random samples from iterable provided
and list of its indexes
'''
   l = len(iterable)
   k = set()
   while n>len(k):
      k.add(randint(0, l-1))

   l = []
   wk = deepcopy(iterable)
   for i in list(k):
      l.append(wk[i])
   return list(k),l


def negativeSampling(sampling,iterable):
   '''subsamples negative of selection list
requires:
   [iterable: str(), tuple(), list()]
   [sampling: tuple(), list()] indexes[int()] *not* to retrieve
ensures:
   collection without indexes present in sampling
'''
   j = deepcopy(list(iterable))
   s = sorted(list(set(sampling)), reverse=True)
   for i in s:
      j.pop(i)
   return j


if __name__=='__main__':
   print('test')
   k = 'asdfg'
   print(WithReposition(4,k))
   print(WithoutReposition(4,k))
   print(negativeSampling((1,4,2,3,9),"0123456789"))


