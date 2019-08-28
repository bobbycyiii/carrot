import regina
import fault

def isCommonAxisCommutator(expr):
    l = expr.terms()
    if len(l) != 4 or \
       l[0] != l[2].inverse() or \
       l[1] != l[3].inverse():
        return False
    else:
        return True

def isCommonAxisEquation(expr):
    l = expr.terms()
    return len(l) == 2

def isCommonAxisPresentation(G):
    adjSet = {}
    for g in range(G.countGenerators()):
        adjSet[g] = set()
    for i in range(G.countRelations()):
        r = G.relation(i)
        a = r.generator(0)
        b = r.generator(1)
        if isCommonAxisCommutator(r) or isCommonAxisEquation(r):
            adjSet[a].add(b)
            adjSet[b].add(a)
    comp = set([0,])
    boundary = adjSet[0]
    while len(boundary) > 0:
        newboundary = set()
        for v in boundary:
            comp.add(v)
        for v in boundary:
            for w in adjSet[v]:
                if not w in comp:
                    newboundary.add(w)
        boundary = newboundary
    return len(comp) == G.countGenerators()

def addDipyramid(T,n):
    N = T.size()
    for i in range(n):
        T.newTetrahedron()
    for i in range(n):
        me = T.tetrahedron(N+i)
        you = T.tetrahedron(N+((i+1)%n))
        me.join(2,you,regina.Perm4(2,3))

def makeDipyramids(sig):
    T = regina.Triangulation3()
    for i in range(len(sig)):
        ni = sig[i]
        addDipyramid(T,ni)
    return T
    
def milleyWhichFace(label,n):
    if label < n:
        tet, face = label, 0
    else:
        tet, face = label-n, 1
    return (tet, face)
    
def milleyGluing(gluing):
    (sig, cayley) = gluing
    n = sum(sig)
    assert len(cayley) == 2*n
    M = makeDipyramids(sig)
    for i in range(2*n):
        (me,myFace) = milleyWhichFace(i,n)
        (you,yourFace) = milleyWhichFace(cayley[i],n)
        if myFace == yourFace:
            sigma = regina.Perm4(2,3)
        else:
            sigma = regina.Perm4(0,1)
        M.tetrahedron(me).join(myFace,M.tetrahedron(you),sigma)
    return M

def milleyGluingToString(gluing):
    outStr = "("
    for i in gluing[0]:
        outStr = outStr + str(i)
    outStr = outStr + "):"
    for i in gluing[1]:
        outStr = outStr + " " + str(i)
    return outStr

def getSig(s):
    (psig,paren,rest) = s.partition(')')
    sig = psig[1:]
    return map(lambda ch:int(ch), sig)
    
def getTuple(s):
    (gluing,paren,rest) = s.rpartition('(')
    (sig,colon,numbers) = gluing.partition(':')
    return map(lambda s: int(s), numbers.split(' ')[1:-1])

def badGroup(pres):
    recog = pres.recogniseGroup()
    if recog != '':
        if recog[0] != 'Z':
            return True
        elif recog == "Z":
            return True
        elif recog[1] != "~":
            return True
    elif G.countRelations() == 0:
        return True
    elif G.identifyAbelian():
        return True
    elif isCommonAxisPresentation(G):
        return True
    else:
        return False

import sys
if __name__ == "__main__":
    filename = sys.argv[1]
    f = open(filename, 'r')
    x = f.readline()
    isoSigs = {}
    hypNames = {}
    nonhyperbolics = 0
    badGroups = 0
    faulty = 0
    while x != '':
        if x[0] == '<':
            x = f.readline()
            continue
        gluing = (getSig(x),getTuple(x))
        M = milleyGluing(gluing)
        M.intelligentSimplify()
        G = M.fundamentalGroup()
        # if isoSig in isoSigs:
        #    pass
        if badGroup(G):
            # The following assertion checks that we get the same result as Milley did for this gluing.
            assert not "geometric" in x
            print milleyGluingToString(gluing) + ": not hyperbolic: common axis presentation"
            badGroups = badGroups + 1
            nonhyperbolics = nonhyperbolics + 1
        elif M.hasStrictAngleStructure():
            # The following assertion checks that we get the same result as Milley did for this gluing.
            assert "geometric" in x
            h = regina.Census.lookup(M)
            j = 0
            while h.empty():
                M.simplifyExhaustive(j)
                h = regina.Census.lookup(M)
                j = j + 1
            hit = h.first().name().split(' ')[0]
            print milleyGluingToString(gluing) + ": is hyperbolic: census: " + hit
            if not hit in hypNames:
                hypNames[hit] = []
            hypNames[hit].append((gluing,M.isoSig()))
        else:
            if fault.isFaultless(M):
                # This never happens for Mom-4s.
                print milleyGluingToString(gluing) + ": is hyperbolic: faultless"
            else:
                # The following assertion checks that we get the same result as Milley did for this gluing.
                assert not "geometric" in x
                print milleyGluingToString(gluing) + ": not hyperbolic: has fault: " + M.isoSig()
                faulty = faulty + 1
                nonhyperbolics = nonhyperbolics + 1
        x = f.readline()
        continue
    hN = list(hypNames)
    hN.sort()
    for hit in hN:
        outStr = hit + " |"
        for p in hypNames[hit]:
            (gluing,isoSig) = p
            outStr = outStr + " " + milleyGluingToString(gluing) + "; " + isoSig + " |"
        print outStr
    print "Results for {0}:".format(filename)
    print "There are {0} hyperbolic Mom-4 manifolds up to homeomorphism.".format(len(hN))
    print "There were {0} nonhyperbolic gluings.".format(nonhyperbolics)
    print "There were {0} common axis presentations and {1} other manifolds with faults.".format(badGroups, faulty)
