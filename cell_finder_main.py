## TODO
# 3/4/20
# archive previous 10x ref images -- IMPLEMENTED, TEST!
# 2/20/20
# batch close previews? or toggle them on/off more easily -- DONE
# 2/13/20
# add a,b,c suffixes to stage names e.g. `A02c` -- DONE
# top of STG file should have num of positions -- NOT DONE, JUST SET TO 999
# label preview image with cell index -- DONE
# add off-center delta adjustment -- DONE


## WOULD BE NICE
# do not image adjacent cells

import os.path, h5py, cv2, csv, re, string
import numpy as np
import matplotlib.pyplot as plt
from importlib import reload
import shutil
# %matplotlib inline

""" USER-SET PARAMETERS"""
# parameters
DEBUGGING = 0 # set to 1 if running debug version (ilya laptop)
showImages = 0 # 1 to show images
moveToArchive = 1 # after processing, move 10X tif files and ilastik h5 files to an archive directory
minArea = 80 # minimum area to be considered a cell (in pixels)
maxArea = 250 # max area to be considered a cell (in pixels)
distToEdgeThresh = 50 # cell must be at least distToEdgeThresh (in pixels) away from FOV edge
maxNumCells = 9 # max number of cells to output
pxToStepsX, pxToStepsY = (1.76, 1.76) # steps x movements = pixels * pxToSteps
imgSize = 512 # edge size of 10x image
centerAdjustmentX, centerAdjustmentY = (52, -34) # correct for offset in steps between 40x and 10x image centers
keyName = 'exported_data'

if DEBUGGING:
    # ILYA LOCAL PATHS
    defaultSTGFileDir = r'sample.STG' # path to default STG file
    outputSTGFileDir = r'sampleOut.STG' # path to STG file that will be saved
    h5FullDir = r'Z:\ilya\projects\cell_finding\10xRefImages_P5_NoCrop' # path to ref images folder (AutoFocusRef1 images)
    archiveDir = r'Z:\ilya\projects\cell_finding\img_archive'
else:
    # RIG 1 PATHS
    defaultSTGFileDir = r'D:\\Neuronal_Culture_DataMM\\refs\\sampleIn.STG' # path to default STG file
    outputSTGFileDir = r'D:\\Neuronal_Culture_DataMM\\refs\\sampleOut.STG' # path to STG file that will be saved
    h5FullDir = r'D:\Neuronal_Culture_DataMM\\refs\\temp' # path to ref images folder (AutoFocusRef1 images)
    archiveDir = r'D:\Neuronal_Culture_DataMM\\refs\\temp\\archive'


""" /USER-SET PARAMETERS"""

def runIlastik(tifFile):
    ilastik_bat_str = r'"C:\\Program Files\\ilastik-1.3.3post1\\run-ilastik.bat"' if DEBUGGING else r'"C:\\Program Files\\ilastik-1.3.3post2\\run-ilastik.bat"'
    
    ilastik_params = r'--headless --project=voltron_cell_finding.ilp --export_source="Simple Segmentation"'
    parentFolder = os.path.dirname(tifFile)
    ilastik_output_params = r'--output_filename_format="' + os.path.join(parentFolder, r'{nickname}_segmentation.h5"')

    ilastik_cmd = ilastik_bat_str + ' ' + ilastik_params + ' ' + ilastik_output_params + ' ' + '"' + tifFile + '"'
    print(ilastik_cmd)
    # response = os.system('dir')
    response = os.system(r'"' + ilastik_cmd + r'"') # launch headless ilastik (remember quotes!)
    print('Returned ' + str(response))

def plotContour(c, h5d):
    """plots contour outlines"""
    ctrImg = np.ones(np.shape(h5d))
    cv2.drawContours(ctrImg, c, -1, (0, 255, 0), 1)
    cv2.imshow("test", ctrImg)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

'''
plots ilastik-generated mask, cell locations (white circle), and cell index numbers
'''
def plotCellLocs(img, x, y):
    fig, ax = plt.subplots(1)
    plt.imshow(img)
    xToPlot, yToPlot = (np.asarray(x, dtype=np.float32), np.asarray(y, dtype=np.float32))
    
    nCells = len(x)
    
    plt.plot(xToPlot, yToPlot, 'wo')
    for cellNum in range(nCells): 
        print(str(cellNum) + ' ' + str(xToPlot[cellNum]) + ' ' + str(yToPlot[cellNum]))
        ax.annotate(str(cellNum), (xToPlot[cellNum]+10, yToPlot[cellNum]+10), color='white')
    ax.set_aspect('equal', 'box')
    ax.set(xlim=(0, 512), ylim=(0, 512))
    plt.gca().invert_yaxis()
    plt.show()


'''
loads ilastik-parsed h5 file, performs blob detection to weed out bad cells, prunes cell list to only good cells
returns XY array of cell centers in pixels
'''
def getCellCoords(h5Filename):
    h5f = h5py.File(h5Filename, 'r')
    h5data = np.abs( np.squeeze(h5f[keyName]) - 2) # binarize to [0,1]
    # print(h5data.shape)
    

    xSize = h5data.shape[0]
    ySize = h5data.shape[1]
    
    # find cell contours
    contours, h = cv2.findContours(h5data, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    
    i=0
    coordsX = np.zeros(len(contours))
    coordsY = np.zeros(len(contours))
    areas = np.zeros(len(contours))
    for c in contours:
        M = cv2.moments(contours[i])
        if M["m00"] != 0: 
            coordsX[i] = int(M["m10"] / M["m00"])
            coordsY[i] = int(M["m01"] / M["m00"])
            areas[i] = cv2.contourArea(c)
        else: # segmentation error
            coordsX[i] = 0
            coordsY[i] = 0
            areas[i] = 0

        i+=1
    
    print("# found cells: " + str(len(contours)))
    # plotContour(contours, h5data)
    
    # remove small spots, cells found on the edge of FOV, etc

    # true to remove
    toRemove = (coordsX < distToEdgeThresh) | (coordsY < distToEdgeThresh) | \
        (xSize - coordsX < distToEdgeThresh) | (ySize - coordsY < distToEdgeThresh) | (areas < minArea) | (areas > maxArea)
    good_Contours = [d for (d, remove) in zip(contours, toRemove) if not remove]
    good_coordsX = [d for (d, remove) in zip(coordsX, toRemove) if not remove]
    good_coordsY = [d for (d, remove) in zip(coordsY, toRemove) if not remove]
    good_areas = [d for (d, remove) in zip(areas, toRemove) if not remove]

    print("# good cells: " + str(len(good_Contours)))

    # sort cells best > worst according to area (largest area = best)
    sorted_good_areas = sorted(good_areas, reverse=1)
    # sorted_good_Contours = [x for _,x in sorted(zip(good_areas, good_Contours), key=lambda p: p[0], reverse=1)]
    
    sorted_good_coordsXY = [(x,y) for _,x,y in sorted(zip(good_areas, good_coordsX, good_coordsY), reverse=1)]
    sorted_good_coordsX  = [j[0] for j in sorted_good_coordsXY]
    sorted_good_coordsY  = [j[1] for j in sorted_good_coordsXY]
    
    # truncate list to maxNumCells
    if len(sorted_good_coordsX) > maxNumCells:
        sorted_good_coordsX = sorted_good_coordsX[0:maxNumCells]
        sorted_good_coordsY = sorted_good_coordsY[0:maxNumCells]

    if showImages:
        plotCellLocs(h5data, sorted_good_coordsX, sorted_good_coordsY)

    print('XY coordinates (pixels):')
    print(np.vstack((np.asarray(sorted_good_coordsX), np.asarray(sorted_good_coordsY))))
    return np.vstack((np.asarray(sorted_good_coordsX), np.asarray(sorted_good_coordsY)))

'''
modifies well titile like "96Well02-A02_500dot456" to "96Well02-A02c_500dot456" (adds letter after well position)
'''
def modifyWellTitle(title, letter):
    wellString = re.split('_', wellTitle)
    return wellString[0] + letter + '_' + wellString[1]
    
# START MAIN SCRIPT

allTifs = [os.path.join(h5FullDir,f) for f in os.listdir(h5FullDir) if f.endswith('.tif')]

if DEBUGGING:
    # FOR DEBUGGING ONLY: fewer tifs
    allTifs = allTifs[:10] 
    #####

print(str(len(allTifs)) + ' TIF files')

if not DEBUGGING:
    # FOR DEBUGGING ONLY: do not run ilastik
    for tif in allTifs: runIlastik(tif)

imgCenter = np.round(imgSize / 2)
letters = string.ascii_lowercase[:27]
# load and modify STG file
# assume STG file contains XYZ positions of center FOV in each well (focused in 40x)
# generate new STG file with scaled positions of found cells
nHeaderLines = 3

allH5files = [os.path.join(h5FullDir,f) for f in os.listdir(h5FullDir) if f.endswith('.h5')]
print(str(len(allH5files)) + ' h5 files')

with open(defaultSTGFileDir) as stgInFile, open(outputSTGFileDir, 'w') as stgOutFile:
    i=0
    
    for row in stgInFile:
        # print(row)
        
        # copy header lines
        if i < nHeaderLines:
            stgOutFile.writelines(row) 
        elif i == nHeaderLines: # this will write the `num FOVs` line. Should actually compute the num of FOVs but may not be necessary
            stgOutFile.writelines('999\n')
        else:
            
            # load corresponding h5, process, get cell indices
            r = row.split(',', 3) # [name, X, Y, rest of string]
            wellTitle = r[0]
            wellID = re.split('-|_', wellTitle)[1] # e.g. 'A05'
            
            # find corresponding h5 file
            match = [s for s in allH5files if wellID in s]
            if len(match) > 1 or len(match) == 0:
                raise ValueError(wellID + ': mistake in finding H5 files')
            
            XYCoords = getCellCoords(match[0])
            Xcoords = XYCoords[0,]
            Ycoords = XYCoords[1,]
            
            nFoundCells = len(Xcoords)
            
            # relative move in steps
            # formula for move:
            #    well center + (desired loc - screen center) * pxToSteps - stepAdjustment
            sorted_good_coordsX_steps_abs = int(r[1]) + (np.array(Xcoords) - imgCenter) * pxToStepsX - centerAdjustmentX
            sorted_good_coordsY_steps_abs = int(r[2]) + (np.array(Ycoords) - imgCenter) * pxToStepsY - centerAdjustmentY
            
            # add row for each found cell
            for j in range(nFoundCells):
                # write STG line e.g.: "96Well02-A02a_500dot456", -16788.28, -5936.32, 8577.79, 0, 8577.79, FALSE, -9999, TRUE, TRUE, 0, -1
                writeString = modifyWellTitle(wellTitle, letters[j]) + ', ' + str(sorted_good_coordsX_steps_abs[j]) + ', ' + str(sorted_good_coordsY_steps_abs[j]) + ',' + r[3]
                stgOutFile.writelines(writeString)
            
            
        i+=1


if moveToArchive:
    print('Archiving tif and h5 segmentation files... ')
    for f in allTifs:
        shutil.move(f, os.path.join(archiveDir, os.path.split(f)[1]))
    for f in allH5files:
        shutil.move(f, os.path.join(archiveDir, os.path.split(f)[1]))
    print('done\n')

print('>>>>>>>>>>>>>>>>>>>>>>SUCCESS!!!!!!!!!!!!<<<<<<<<<<<<<<<<<<<<<<<<<<')