import numpy as np

# TODO Add ability to name nodes/edges
# TODO Add error checking
# TODO Add Network Completion

# Rows correspond to vertices
# Columns correspond to edges


class Ubergraph:

    # Construct from list of edges
    def __init__(self, *incidence):
        self.matrix = np.column_stack(incidence)
        self.order = np.shape(self.matrix)[0]
        self.size = np.shape(self.matrix)[1]
        self.nodeNames = []
        self.edgeNames = []
        for i in range(self.order):
            self.nodeNames.append("v_" + str(i))
        for i in range(self.size):
            self.edgeNames.append("e_" + str(i))

    def getOrder(self):
        return self.order

    def getSize(self):
        return self.size

    def getNodeName(self, id):
        return self.nodeNames[id]

    def getEdgeName(self, id):
        return self.edgeNames[id]

    def getNodeID(self, name):
        return self.nodeNames.index(name)

    def getEdgeID(self, name):
        return self.edgeNames.index(name)

    def setNodeName(self, id, newName):
        self.nodeNames[id] = newName

    def setEdgeName(self, id, newName):
        self.edgeNames[id] = newName

    def getNode(self, id):
        return self.matrix[id, :]

    def setNode(self, id, incidence):
        self.matrix[id, :] = incidence

    def getEdge(self, id):
        return self.matrix[:, id]

    def setEdge(self, id, incidence):
        self.matrix[:, id] = incidence

    def addEdge(self, incidence):
        self.matrix = np.concatenate((self.matrix, incidence.T), axis=1)
        self.size = self.size + 1

    def addNode(self, incidence):
        self.matrix = np.concatenate((self.matrix, incidence), axis=0)
        self.order = self.order + 1

    # print matrix
    def print(self):
        print(self.matrix)


if __name__ == "__main__":
    uber = Ubergraph(
        np.array([1, 0, 0, 0]), np.array([0, -2, 0, 0]), np.array([1, 1, 0, 0])
    )
    uber.print()
    print("order: " + str(uber.getOrder()))
    print("size: " + str(uber.getSize()))
    uber.addEdge(np.array([[-1, -1, -1, -1]]))
    uber.print()
    uber.addNode(np.array([[8, 8, 8, 8]]))
    uber.print()
    print("\n")
    print(uber.getEdge(1))
    print("order: " + str(uber.getOrder()))
    print("size: " + str(uber.getSize()))
    print("\n")
    print(uber.getNodeName(2))
    uber.setNodeName(2, "Pebus")
    print(uber.getNodeName(2))
    print(uber.getNodeID("Pebus"))
    print(uber.getEdgeName(1))
