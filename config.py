import os

if os.path.isdir('/TL/euresist_phylodynamics/work/Projects/MSTBasedForests/scripts/'):
    projectPath='/TL/euresist_phylodynamics/work/Projects/MSTBasedForests/'    
    toolPath = '/TL/euresist_phylodynamics/work/Tools/'
    scriptPath = '/TL/euresist_phylodynamics/work/Projects/MSTBasedForests/scripts/'
    tempPath = '/TL/euresist_phylodynamics/work/Projects/MSTBasedForests/temp/'
    debugPath = '/TL/euresist_phylodynamics/work/Projects/MSTBasedForests/debug/'
elif os.path.isdir('/local/home/pk/Projects/MSTBasedForests/scripts/'):    
    projectPath = '/local/home/pk/Projects/MSTBasedForests/'
    scriptPath = '/home/pk/workspace/MSTBasedForests/'
    toolPath = '/usr/local/bin/'
    tempPath = '/local/home/pk/Projects/MSTBasedForests/temp/'
    playgroundPath = '/home/pk/playground/'
    treeShapePath = '/home/pk/Projects/mstBasedPhylogenetics/results/interestingShapes/'
    debugPath = '/local/home/pk/Projects/MSTBasedForests/debug/'
elif os.path.isdir('/project/exaptation/Projects/MSTBasedForests/scripts/'):
    projectPath='/project/exaptation/Projects/MSTBasedForests/'    
    toolPath = '/project/exaptation/Tools/'
    scriptPath = '/project/exaptation/Projects/MSTBasedForests/scripts/'
    tempPath = '/project/exaptation/Projects/MSTBasedForests/temp/'
    debugPath = '/project/exaptation/Projects/MSTBasedForests/debug/'
elif os.path.isdir('/home/kalaghat/exaptation/Projects/MSTBasedForests/scripts/'):
    projectPath='/home/kalaghat/exaptation/Projects/MSTBasedForests/'    
    toolPath = '/home/kalaghat/exaptation/Tools/'
    scriptPath = '/home/kalaghat/exaptation/Projects/MSTBasedForests/scripts/'
    tempPath = '/home/kalaghat/exaptation/Projects/MSTBasedForests/temp/'
    debugPath = '/home/kalaghat/exaptation/Projects/MSTBasedForests/debug/'
elif os.path.isdir('/home/kalaghat/Projects/workspace_scripts/'):
    projectPath='/home/kalaghat/Projects/proj_dir/'    
    toolPath = '/home/kalaghat/Projects/Tools/'
    scriptPath = '/home/kalaghat/Projects/workspace_scripts/'
    tempPath = '/home/kalaghat/Projects/temp/'
    debugPath = '/home/kalaghat/Projects/debug/'
    
 
def GetmxqsubPrefix(groupId=""):
    mxqsubPrefix=""
    if 'exaptation' in os.getcwd():
        mxqsubPrefix = "mxqsub -t 12h -m 2G "
        if groupId != "":
            mxqsubPrefix += str(groupId)+" " 
    return (mxqsubPrefix) 
 
def GetReportedVariableRange(variable):
    if variable=='treeType':
        variableRange = ['balanced','unbalanced','random']
    elif variable=='samplingProportion':
        variableRange = ['0.0','0.25','0.75','1.0']
    elif variable=='numberOfObservedVertices':
        variableRange = ['100','10000','100000','1000000']
    elif variable=='branchLength':
        variableRange = ['0.004','0.064','0.256']
    elif variable=='sequenceLength':
        variableRange = ['250','500','2000','4000']
    elif variable=='contractedEdgeType':
        variableRange = ['terminal','internal','SA']
    elif variable=='scalingFactor':
        variableRange = ['0.1','1','10','100']
    elif variable=='growthRate':
        variableRange = ['1','5']
    return variableRange
