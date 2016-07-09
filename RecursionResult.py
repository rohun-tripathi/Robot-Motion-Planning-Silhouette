class VectorsToLink :

    def __init__(self):
        self.intersectionPointsAlongOtherAxis = []
        self.criticalPointsToLinkToFuture = []
        self.criticalPointsToLinkToPast =[]

    def setCriticalPointsToLinkToPast(self, criticalPointsToLinkToPast):
        self.criticalPointsToLinkToPast = criticalPointsToLinkToPast
        return self

    def getCriticalPointsToLinkToPast(self):
        return self.criticalPointsToLinkToPast

    def setCriticalPointsToLinkToFuture(self, criticalPointsToLinkToFuture):
        self.criticalPointsToLinkToFuture = criticalPointsToLinkToFuture
        return self

    def getCriticalPointsToLinkToFuture(self):
        return self.criticalPointsToLinkToFuture

    def setIntersectionPointsAlongOtherAxisReturnSelf(self, intersectionPointsAlongOtherAxis):
        self.intersectionPointsAlongOtherAxis = intersectionPointsAlongOtherAxis
        return self

    def getIntersectionPointsAlongOtherAxis(self):
        return self.intersectionPointsAlongOtherAxis