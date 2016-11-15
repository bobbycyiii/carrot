import regina
def possiblyHyp(mfld):
  m = mfld
  return m.isValid() \
     and m.isOrientable() \
     and m.getEulerCharManifold() == 0 \
     and m.isConnected()

def separates(surf):
  M = surf.cutAlong()
  return not M.isConnected()

def isAnnulus(surf):
  return surf.hasRealBoundary()            \
    and surf.getEulerCharacteristic() == 0 \
    and surf.isOrientable()

def isNonSeparatingAnnulus(surf):
  return isAnnulus(surf) \
     and not separates(surf)

def isAnnulusFault(surf):
  return isAnnulus(surf) and \
         isFault(surf)

def findNonSeparatingAnnulus(mfld):
  T = mfld
  nsl = regina.NNormalSurfaceList.enumerate
  std = regina.NS_STANDARD
  fnd = regina.NS_FUNDAMENTAL
  l = nsl(T,std,fnd)
  a = None
  n = l.getNumberOfSurfaces()
  for i in range(0,n):
    surf = l.getSurface(i)
    if isNonSeparatingAnnulus(surf):
      a = surf
      break
  return a

def embedded(edge):
  src = edge.getVertex(0)
  snk = edge.getVertex(1)
  return src != snk

def lrmaps(edge):
  embs = edge.getEmbeddings()
  return (embs[0].getVertices(),\
          embs[-1].getVertices())

def lrtets(edge):
  embs = edge.getEmbeddings()
  return (embs[0].getTetrahedron(),\
          embs[-1].getTetrahedron())

def foldable(edge):
  if not edge.isBoundary():
    return False
  (F0,F_1) = lrmaps(edge)
  (D0,D_1) = lrtets(edge)
  lvx = D0.getVertex(F0[2])
  rvx = D_1.getVertex(F_1[3])
  return lvx != rvx

def twoTwo(edge):
  M = edge.getTriangulation()
  T = M.newTetrahedron()
  (F0,F_1) = lrmaps(edge)
  (D0,D_1) = lrtets(edge)
  X = regina.NPerm(2,3)
  S0 = F0 * X
  S_1 = F_1 * X
  T.joinTo(2,D0,S0)
  T.joinTo(3,D_1,S_1)

def foldAlong(edge):
  assert edge.isBoundary()
  (D0,D_1) = lrtets(edge)
  (F0,F_1) = lrmaps(edge)
  X = regina.NPerm(2,3)
  glu = F_1 * X * F0.inverse()
  D0.joinTo(F0[3], D_1, glu)

def firstBoundaryEdge(mfld,pred):
  cpts = mfld.getBoundaryComponents()
  for d in cpts:
    n = d.getNumberOfEdges()
    for i in range(0,n):
      e = d.getEdge(i)
      if pred(e):
        return e
  else:
    return None

def simplifyCusps(finite_mfld):
  M = finite_mfld
  fBE = firstBoundaryEdge
  f = fBE(M,embedded)
  while f != None:
    twoTwo(f)
    g = fBE(M,foldable)
    while g != None:
        foldAlong(g)
        g = fBE(M,foldable)
    f = fBE(M,embedded)

def possiblyT2xI(mfld):
  m = mfld
  dno = m.getNumberOfBoundaryComponents()
  if not (possiblyHyp(m) and dno == 2):
    return False
  cpts = m.getBoundaryComponents()
  for d in cpts:
    x = d.getEulerCharacteristic()
    if x != 0:
      return False
  H1 = m.getHomologyH1()
  if not H1.toString() == '2 Z':
    return False
  H2 = m.getHomologyH2()
  if not H2.isZ():
    return False
  H1R = m.getHomologyH1Rel()
  if not H1R.isZ():
    return False
  return True

def hasTorusBoundary(mfld):
  m = mfld
  if m.isClosed():
    return False
  for v in m.getVertices():
    if not v.getLink() == regina.NVertex.TORUS:
      return False
  return True

def irreducible(regina_mfld):
  M = regina.NTriangulation(regina_mfld)
  M.finiteToIdeal()
  M.intelligentSimplify()
  nsl = regina.NNormalSurfaceList.enumerate
  l = nsl(M, regina.NS_STANDARD,  \
             regina.NS_FUNDAMENTAL)
  n = l.getNumberOfSurfaces()
  for i in range(0,n):
    s = l.getSurface(i)
    x = s.getEulerCharacteristic()
    if x != 2:
      continue
    if isFault(s):
      return False
  else:
    return True

def isT2xI(regina_mfld):
  M = regina.NTriangulation(regina_mfld)
  M.finiteToIdeal()
  M.intelligentSimplify()
  if not possiblyT2xI(M):
    return False
  D = regina.NTriangulation(M)
  simplifyCusps(D)
  T = D.getBoundaryComponent(1)
  n = T.getNumberOfEdges()
  assert n == 3
  for i in range(0,n):
    clone = regina.NTriangulation(D)
    cpt = clone.getBoundaryComponent(1)
    e = cpt.getEdge(i)
    foldAlong(e)
    if not clone.isSolidTorus():
      return False
  else:
    return True

b = lambda m: m.isBall()
cd = lambda m: m.hasCompressingDisc()
d2s1 = lambda m: m.isSolidTorus()
t2i = isT2xI
def isFault(surf):
  s = surf
  x = s.getEulerCharacteristic()
  definitelyNotFault =     \
    x < 0 or               \
    not s.isConnected() or \
    s.isVertexLink()
  if definitelyNotFault:
    return False
  if not s.isOrientable():
    return True
  M1 = s.cutAlong()
  M1.intelligentSimplify()
  if M1.isConnected():
    return True
  assert M1.splitIntoComponents() == 2
  M1.intelligentSimplify()
  M2 = M1.getFirstTreeChild()
  M3 = M2.getNextTreeSibling()
  if s.hasRealBoundary():
    if x == 1:
      # s is a disc
      return b(M2) == b(M3)
    else:
      # s had better be an annulus
      assert x == 0
      return not (b(M2) or b(M3)) \
         and d2s1(M2) == d2s1(M3)
  else:
    # s is closed
    if x == 2:
      # s is a sphere
      return not (b(M2) or b(M3))
    else:
      # s had better be a torus
      assert x == 0
      return not (t2i(M2) or t2i(M3) \
             or   cd(M2)  or cd(M3))

def unhypByNormalSurfaces(mfld):
  dno = mfld.getNumberOfBoundaryComponents()
  assert dno > 0
  if not possiblyHyp(mfld):
    return True
  T = regina.NTriangulation(mfld)
  T.finiteToIdeal()
  T.intelligentSimplify()
  nsl = regina.NNormalSurfaceList.enumerate
  std = regina.NS_STANDARD
  fnd = regina.NS_FUNDAMENTAL
  l = nsl(T,std,fnd)
  n = l.getNumberOfSurfaces()
  for i in range(0,n):
    surf = l.getSurface(i)
    if isFault(surf):
      return True
  TT = regina.NTriangulation(mfld)
  TT.idealToFinite()
  TT.intelligentSimplify()
  if TT.hasCompressingDisc():
    return True
  vtx = regina.NS_VERTEX
  qd  = regina.NS_QUAD
  ll = nsl(TT,qd,vtx)
  nn = ll.getNumberOfSurfaces()
  for ii in range(0,nn):
    surf = ll.getSurface(ii)
    if TT.getNumberOfBoundaryComponents() == 2:
      if isNonSeparatingAnnulus(surf):
        return True
    else:
      if isAnnulusFault(surf):
        return True
  else:
    return False

def isHyp(regina_manifold):
  m = regina.NTriangulation(regina_manifold)
  m.intelligentSimplify()
  if m.hasStrictAngleStructure():
    return True
  else:
    return not unhypByNormalSurfaces(m)

