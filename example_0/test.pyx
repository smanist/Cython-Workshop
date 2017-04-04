import cython
cimport cython

def getAnswer():
    cdef int answer = 42
    print("\nMy answer is {0}!\n".format(answer))
