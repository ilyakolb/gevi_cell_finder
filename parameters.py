""" USER-SET PARAMETERS"""
DEBUGGING = 1 # set to 1 if running debug version (ilya laptop)
showImages = 0 # 1 to show images
moveToArchive = 1 # after processing, move 10X tif files and ilastik h5 files to an archive directory
minArea = 80 # minimum area to be considered a cell (in pixels)
maxArea = 250 # max area to be considered a cell (in pixels)
distToEdgeThresh = 50 # cell must be at least distToEdgeThresh (in pixels) away from FOV edge
maxNumCells = 9 # max number of cells to output
pxToStepsX, pxToStepsY = (1.76, 1.76) # steps x movements = pixels * pxToSteps
imgSize = 512 # edge size of 10x image
centerAdjustmentX, centerAdjustmentY = (52, -34) # correct for offset in steps between 40x and 10x image centers
nHeaderLines = 3 # number of header lines in STG file

# ilastik params
ilastik_bat_str = r'"C:\\Program Files\\ilastik-1.3.3post1\\run-ilastik.bat"' if DEBUGGING else r'"C:\\Program Files\\ilastik-1.3.3post2\\run-ilastik.bat"'
ilastik_project_name = 'voltron_cell_finding.ilp' # name of ilastik project to run
keyName = 'exported_data' # this is set by ilastik

# note: h5 files get saved in the same directory as original tif files
if DEBUGGING:
    # ILYA LOCAL PATHS
    defaultSTGFileDir = r'sample.STG' # path to default STG file
    outputSTGFileDir = r'sampleOut.STG' # path to STG file that will be saved
    h5FullDir = r'testing' # path to ref images folder (AutoFocusRef1 images)
    archiveDir = r'img_archive'
else:
    # RIG 1 PATHS
    defaultSTGFileDir = r'D:\\Neuronal_Culture_DataMM\\refs\\sampleIn.STG' # path to default STG file
    outputSTGFileDir = r'D:\\Neuronal_Culture_DataMM\\refs\\sampleOut.STG' # path to STG file that will be saved
    h5FullDir = r'D:\Neuronal_Culture_DataMM\\refs\\temp' # path to ref images folder (AutoFocusRef1 images)
    archiveDir = r'D:\Neuronal_Culture_DataMM\\refs\\temp\\archive'


""" /USER-SET PARAMETERS"""