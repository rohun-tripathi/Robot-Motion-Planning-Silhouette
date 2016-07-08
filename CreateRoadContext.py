import shared as SHARED

class RoadContext :

    def __init__(self):
        self.traversalAxis = 0
        self.criticalPoints = []
        self.parentVector = [0 for x in range(SHARED.dim)]
        self.originList = []
        self.ellipseList = []

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