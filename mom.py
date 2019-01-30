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

def isCommonAxisPresentation(pres):
    if G.countGenerators() == 2:
        for i in range(G.countRelations()):
            if isCommonAxisCommutator(G.relation(i)):
                return True
            if isCommonAxisEquation(G.relation(i)):
                return True
    return False

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
    
def whichFace(label,n):
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
        (me,myFace) = whichFace(i,n)
        (you,yourFace) = whichFace(cayley[i],n)
        if myFace == yourFace:
            sigma = regina.Perm4(2,3)
        else:
            sigma = regina.Perm4(0,1)
        M.tetrahedron(me).join(myFace,M.tetrahedron(you),sigma)
    return M

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
    while x != '':
        if x[0] == '<':
            x = f.readline()
            continue
        gluing = (getSig(x),getTuple(x))
        M = milleyGluing(gluing)
        M.intelligentSimplify()
        isoSig = M.isoSig()
        G = M.fundamentalGroup()
        if isoSig in isoSigs:
            pass
        elif badGroup(G):
            isoSigs[isoSig] = str(gluing) + ": not hyp: bad group"
        elif M.hasStrictAngleStructure():
            h = regina.Census.lookup(M)
            x = 0
            while h.empty():
                M.simplifyExhaustive(x)
                h = regina.Census.lookup(M)
                x = x + 1
            else:
                hit = h.first().name().split(' ')[0]
                sigOut = str(gluing) + ": hyperbolic: census"
                isoSigs[isoSig] = sigOut + ": " + hit
                if not hit in hypNames:
                    hypNames[hit] = []
                hypNames[hit].append((gluing,isoSig))
        elif fault.isFaultless(M):
            isoSigs[isoSig] = str(gluing) + ": not hyp: has fault"
        else:
            isoSigs[isoSig] = str(gluing) + ": not hyp: faultless"
        print isoSigs[isoSig]
        x = f.readline()
        continue
    hN = list(hypNames)
    hN.sort()
    for hit in hN:
        outStr = hit + "|"
        for p in hypNames[hit]:
            (gluing,isoSig) = p
            (sig,cayley) = gluing
            outStr = outStr + " ("
            for i in sig:
                outStr = outStr + str(i)
            outStr = outStr + "):"
            for i in cayley:
                outStr = outStr + " " + str(i)
            outStr = outStr + "; " + isoSig + " |"
        print outStr
