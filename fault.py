# Python code accompanying hyperbolicity test
import regina

def isFaultless(T):
    # Validity follows from orientability.
    # assert isMaterial(T) 
    # assert T.isOrientable()
    assert T.hasRealBoundary()
    # If T is ideal, this routine can segfault!
    
    nsl = regina.NormalSurfaces.enumerate
    std = regina.NS_QUAD
    fnd = regina.NS_VERTEX
    l = nsl(T,std,fnd)
    return not hasFault(l)

def hasFault(nsl):
    T = regina.Triangulation3(nsl.triangulation())
    def has(pred):
        N = nsl.size()
        for i in range(N):
            x = nsl.surface(i)
            if pred(x):
                return True
        else:
            return False
    if has(smallNonOrientable) \
       or has(essentialS2) \
       or T.hasCompressingDisc() \
       or has(essentialT2) \
       or has(seifertA2):
        return True
    else:
        return False
    
def smallNonOrientable(surf):
    if surf.eulerChar() >= 0 \
       and not surf.isOrientable():
        return True
    else:
        return False

def essentialS2(surf):
    if surf.eulerChar() != 2:
        return False
    X = chop(surf)
    if len(X) == 1:
        return True
    (LL, RR) = X
    L = regina.Triangulation3(LL)
    R = regina.Triangulation3(RR)
    if not (L.isBall() or R.isBall()):
        return True
    else:
        return False

def chop(surf):
    X = surf.cutAlong()
    X.intelligentSimplify()
    if X.isConnected():
        return (X,)
    X.splitIntoComponents()
    left = X.firstChild()
    right = left.nextSibling()
    return (left.isoSig(),right.isoSig())

def essentialT2(surf):
    if surf.eulerChar() != 0 \
       or not surf.isOrientable() \
       or surf.hasRealBoundary():
        return False
    X = chop(surf)
    if len(X) == 1:
        return True
    (LL,RR) = X
    L = regina.Triangulation3(LL)
    R = regina.Triangulation3(RR)
    t2i = lambda y: isT2xI(y)
    if not (L.hasCompressingDisc() or R.hasCompressingDisc() \
            or t2i(L) or t2i(R)):
        return True
    else:
        return False

def seifertA2(surf):
    if surf.eulerChar() != 0 \
       or not surf.isOrientable() \
       or not surf.hasRealBoundary():
        return False
    X = chop(surf)
    if len(X) == 1:
        return True
    (LL,RR) = X
    L = regina.Triangulation3(LL)
    R = regina.Triangulation3(RR)
    d2s1 = lambda x: x.isSolidTorus()
    if not (L.isBall() or R.isBall()) \
       and d2s1(L) == d2s1(R):
        return True
    else:
        return False

def clone(T):
    X = regina.Triangulation3(T)
    return X

def isT2xI(T):
    if not isHomologyT2xI(T):
        return False
    X = clone(T)
    X.intelligentSimplify()
    simplifyBoundary(X)
    k = X.boundaryComponent(0)
    for e in k.faces(1):
        i = e.index()
        Xe = clone(X)
        Xe.closeBook(Xe.face(1,i))
        Xe.intelligentSimplify()
        if not Xe.isSolidTorus():
            return False
    else:
        return True

def isHomologyT2xI(T):
    X = clone(T)
    X.idealToFinite()
    X.intelligentSimplify()
    if not X.countBoundaryComponents() == 2:
        return False
    bg = X.homologyBdry()
    if not bg.rank() == 4:
        return False
    h1 = X.homologyH1()
    if not h1.rank() == 2:
        return False
    h1b = X.homologyRel()
    if not h1b.isTrivial():
        return False
    else:
        return True

def simplifyBoundary(T):
    ebe = embeddedBoundaryEdge
    cbe = coembeddedBoundaryEdge
    e = ebe(T)
    while e != None:
        T.layerOn(e)
        f = cbe(T)
        while f != None:
            T.closeBook(f)
            f = cbe(T)
        e = ebe(T)

def embeddedBoundaryEdge(T):
    for k in T.boundaryComponents():
        for e in k.faces(1):
            if e.face(0,0).index() != e.face(0,1).index():
                return e
    else:
        return None

def coembeddedBoundaryEdge(T):
    for k in T.boundaryComponents():
        for f in k.faces(1):
            if T.closeBook(f,perform=False):
                return f
    else:
        return None

