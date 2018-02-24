import math

#coordArray is an array of - [coord1,coord2,...]
#cAlphaPairs is an array of ordered pairs - [[0,1],[1,2],...] - where each ordered pair represents the cAlpha-cAlpha distance that we are creating the pdf for. 
#   For the sake of simplicty, in each pair [a,b] let a<b. This make sure that our pairs don't repeat.

def objectiveFunction(coordArray,cAlphaPairs,templateValues):
    objectiveFunctionSum = 0
    for pair in cAlphaPairs:
        coord1 = coordArray[pair[0]]
        coord2 = coordArray[pair[1]]
        mean = templateValues[pair[0]][pair[1]]
        std = .1
        objectiveFunctionSum = objectiveFunctionSum + math.log(pdf(coord1,coord2,mean,std))
    
    return -objectiveFunctionSum

#coord1,coord2 are [x,y,z] array
def pdf(coord1,coord2,mean,std):
	return 1/math.sqrt(2*math.pi*std**2)*math.exp(-(dist(coord1,coord2)-mean)**2/(2*std**2))

def dist(coord1,coord2):
    return math.sqrt((coord1[0]-coord2[0])**2+(coord1[1]-coord2[1])**2+(coord1[2]-coord2[2])**2)
    
#Gives the gradient of the objective function for the x for the ith carbonAlpha    
def gradientX(i,coordArray,cAlphaPairs,templateValues):
    return gradientOfIthCoordinate(i,0,coordArray,cAlphaPairs,templateValues)

def gradientY(i,coordArray,cAlphaPairs,templateValues):
    return gradientOfIthCoordinate(i,1,coordArray,cAlphaPairs,templateValues)

def gradientZ(i,coordArray,cAlphaPairs,templateValues):
    return gradientOfIthCoordinate(i,2,coordArray,cAlphaPairs,templateValues)

def gradientOfIthCoordinate(i,XYOrZ,coordArray,cAlphaPairs,templateValues):
    gradientSum=0;
    for pair in cAlphaPairs:
        #Only take the derivative if the objective changes with respect to our ith coordinate
        if (pair[0]==i or pair[1]==i):
            coord1 = coordArray[pair[0]]
            coord2 = coordArray[pair[1]]
            mean = templateValues[pair[0]][pair[1]]
            std = .1
            k1 = 1/math.sqrt(2*math.pi*std**2)
            k2 = -1/(2*std**2)
            lastTerm = 1
            if (pair[1]==i):
                lastTerm = -1
            distDerivative = 1/2*1/(dist(coord1,coord2))*2(coord1[XYOrZ]-coord2[XYOrZ])*lastTerm
            gradientOfPdf = 1/(pdf(coord1,coord2,mean,std))*k1*k2*math.exp(k2*(dist(coord1,coord2)-mean)**2)*2*(dist(coord1,coord2)-mean)*distDerivative
            gradientSum = gradientSum + gradientOfPdf
            
    
    
    