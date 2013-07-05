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

# take tableau return tabmac compatable string
def printTab(T):
    l = []
    length = len(T)
    for row in T:
        newRow = []
        for entry in row:
            newRow.append(str(entry))
        if len(l) == len(T)-1:
            l.append(''.join(newRow))
        else:
            l.append(''.join(newRow)+',')
    return '\\scriptsize\\young('+''.join(l)+')'


def bigTable(n):
    ST = StandardTableaux(n)
    pairs = []
    for p in range(1,n+1):
        for q in range(p+1,n+1):
            pairs.append([p,q])
    print('\\begin{center}')
    cols = ['|l|']
    for i in pairs:
        cols.append('l|')
    print('\\begin{longtable}{'+"".join(cols)+'}')
    print('\\hline')
    toprow = []
    for pq in pairs:
        [p,q] = pq
        [s,t] = [str(p),str(q)]
        toprow.append('$s_{%s %s}$'%(s,t))
    print('& '+' & '.join(toprow)+' \\\\')
    print('\\hline')

    for T in ST:
        printT = printTab(T)
        tableRow = []
        tableRow.append(printT)
        for pq in pairs:
            tableRow.append(printTab(spqST(pq,T)))
        print(' & '.join(tableRow)+' \\\\')

    print('\\hline')
    print('\\end{longtable}')
    print('\\end{center}')

# given an spq this defines a permutation action on the set of standard tableau,
# we can calculate what the image of each spq is. There is probably a smarter way
# to do this.
def getPerm(pq,n):
    ST = StandardTableaux(n)
    # ST is not exactly a list so index doesnt work with it. We create a list:
    STlist = []
    for T in ST:
        STlist.append(T)
    m = len(ST)
    perm = []
    for i in range(0,m):
        perm.append(STlist.index(spqST(pq,ST[i]))+1)
    return Permutation(perm)

# Get permutation group of whole J_n in S_m
def wholePerm(n):
    pairs = []
    for p in range(1,n+1):
        for q in range(p+2,n+1):
            pairs.append([p,q])
    spqImages = []
    for pq in pairs:
        spqImages.append(getPerm(pq,n))
    return PermutationGroup(spqImages)

# Get permutation group of only length 2 elements
def len2Perm(n):
    pairs = []
    for p in range(1,n-1):
        pairs.append([p,p+2])
    spqImages = []
    for pq in pairs:
        spqImages.append(getPerm(pq,n))
    return PermutationGroup(spqImages)

def testEq(n):
    N = len2Perm(n)
    G = wholePerm(n)
    print(n)
    if N == G:
        print('they are equal!')
    else:
        print('not equal')
        if N.is_normal(G):
            print('its normal!')
            print(G.quotient(N))


# above we considered J_n acting on all standard tableaux of any shape
# lets look at what happens when we fix a shape
def getPermShape(pq,part):
    ST = StandardTableaux(part)
    # ST is not exactly a list so index doesnt work with it. We create a list:
    STlist = []
    for T in ST:
        STlist.append(T)
    m = len(ST)
    perm = []
    for i in range(0,m):
        perm.append(STlist.index(spqST(pq,ST[i]))+1)
    return Permutation(perm)

# Get permutation group of whole J_n in S_m
def wholePermShape(part):
    n = sum(part)
    pairs = []
    for p in range(1,n+1):
        for q in range(p+2,n+1):
            pairs.append([p,q])
    spqImages = []
    for pq in pairs:
        spqImages.append(getPermShape(pq,part))
    return PermutationGroup(spqImages)

# Get permutation group of only length 2 elements
def len2PermShape(part):
    n = sum(part)
    pairs = []
    for p in range(1,n-1):
        pairs.append([p,p+2])
    spqImages = []
    for pq in pairs:
        spqImages.append(getPermShape(pq,part))
    return PermutationGroup(spqImages)

def testEqShape(part):
    N = len2PermShape(part)
    G = wholePermShape(part)
    if N == G:
        print('they are equal!')
    else:
        print('not equal')
        if N.is_normal(G):
            print('its normal!')
            print(G.quotient(N))

def testall(n):
    parts = Partitions(n)
    for part in parts:
        print('for the partition %s'%(part))
        testEqShape(part)

