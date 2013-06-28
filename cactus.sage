# Code to implement action of cactus group on standard tableaux via action on the crystal B^\otimes n.

from sage.combinat.permutation import robinson_schensted_inverse

# function to create list of lists of zeros for given shape

def createZeroTableau(l):
    master = []
    for n in l:
        row = []
        for i in range(0,n):
            row.append(0)
        master.append(row)
    return master

# function to pop out first entry of tableau, apply jeu de taquin (schutzenberger move and return new tableau old entry and coords of missing box

def SchutzMove(T):
    topRow = T[0]
    restOfT = Tableau(T[1:])
    newT = restOfT.slide_multiply(Tableau([topRow[1:]]))
    coords = []
    if len(newT) < len(T):
        coords = [len(T)-1,0]
    else:
        for r in range(0,len(T)):
            if len(newT[r]) < len(T[r]):
                coords = [r,len(T[r])-1]
    return [newT,T[0][0],coords]


# Sch\"utzenberger involution on tableaux (defined as swapping number around as usual and then adding to all the negatives (where n is the rank of gl_n)
def SchutzenbergerTableaux(T,n):
    ST = createZeroTableau(T.shape())
    while len(T[0])>0:
        [T,i,coords] = SchutzMove(T)
        ST[coords[0]][coords[1]] = n-i+1
    return Tableau(ST)

# Sch\"utzenberger involution on strings 1...n in the tensor product B^n. Input needs to be a list!
def SchutzenbergerCrystal(s,n):
    [P,Q] = Permutation(s[::-1], check_input = False).robinson_schensted()
    SP = SchutzenbergerTableaux(P,n)
    return robinson_schensted_inverse(SP,Q)[::-1]

# now we define the sigmas (cactus commutors)
def cactus(a,b,n):
    Sa = SchutzenbergerCrystal(a,n)
    Sb = SchutzenbergerCrystal(b,n)
    return SchutzenbergerCrystal(Sb+Sa,n)

def cactuspq(pq,s,n):
    p = pq[0]
    q = pq[1]
    if p >= q:
        print('hey! p should be strictly less than q!')
    else:
        firstpart = s[:p-1]
        lastpart = s[q:]
        pth = s[p-1:p]
        middle = s[p:q]
        newMiddle = cactus(pth,middle,n)
        return firstpart+newMiddle+lastpart

def spq(pq,s,n):
    p = pq[0]
    q = pq[1]
    if p >= q:
        print('hey! p should be strictly less than q!')
    else:
        if p == q-1:
            return cactuspq(pq,s,n)
        else:
            return cactuspq(pq,spq([p+1,q],s,n),n)


# now we can calculate the action of s_pq \in J_n on a standard tableaux with n boxes:
def spqST(pq,T):
    s = robinson_schensted_inverse(T,T)
    n = T.size()
    Ss = spq(pq,s,n)
    ST = Permutation(Ss).robinson_schensted()
    return ST[1]

