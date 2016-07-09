import shared as SHARED

class RoadContext :

    def __init__(self):
        self.traversalAxis = 0
        self.criticalPoints = []
        self.parentVector = [0 for x in range(SHARED.dim)]
        self.originList = []
        self.ellipseList = []
        self.sliceEllipseStateList = [0 for x in range(0, len(self.ellipseList))]
        self.sliceVector = []
        self.sliceVectorList = []

    def setSliceVectorListReturnSelf(self, sliceVectorList):
        self.sliceVectorList = sliceVectorList
        return self

    def getSliceVectorList(self):
        return self.sliceVectorList

    def setSliceVectorReturnSelf(self, sliceVector):
        self.sliceVector = sliceVector
        return self

    def getSliceVector(self):
        return self.sliceVector

    def setTraversalAxisReturnSelf(self, traversalAxis):
        self.traversalAxis = traversalAxis
        return self

    def getTraversalAxis(self):
        return self.traversalAxis

    def setCriticalPointsReturnSelf(self, criticalPoints):
        self.criticalPoints = criticalPoints
        return self

    def getCriticalPoints(self):
        return self.criticalPoints

    def setParentVectorReturnSelf(self, parentVector):
        self.parentVector = parentVector
        return self

    def getParentVector(self):
        return self.parentVector

    def setOriginListReturnSelf(self, originList):
        self.originList = originList
        return self

    def getOriginList(self):
        return self.originList

    def setEllipseListReturnSelf(self, ellipseList):
        self.ellipseList = ellipseList
        return self

    def getEllipseList(self):
        return self.ellipseList

    def setSliceEllipseStateListReturnSelf(self, sliceEllipseStateList):
        self.sliceEllipseStateList = sliceEllipseStateList
        return self

    def getSliceEllipseStateList(self):
        return self.sliceEllipseStateList

