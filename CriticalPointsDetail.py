class CriticalPointPair:
    def __init__(self):
        self.firstCP = []
        self.secondCP = []
        self.validityOfFirst = False
        self.validityOfSecond = False

    def setFirstCPReturnSelf(self, firstCP):
        self.firstCP = firstCP
        return self

    def getFirstCP(self):
        return self.firstCP

    def setSecondCPReturnSelf(self, secondCP):
        self.secondCP = secondCP
        return self

    def getSecondCP(self):
        return self.secondCP

    def setValidityOfFirst(self, validityOfFirst):
        self.validityOfFirst = validityOfFirst
        return self

    def getValidityOfFirst(self):
        return self.validityOfFirst

    def setValidityOfSecond(self, validityOfSecond):
        self.validityOfSecond = validityOfSecond
        return self

    def getValidityOfSecond(self):
        return self.validityOfSecond
