pro cLineGapFiller_v2_trimmer

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
;cLineGapFiller.pro
; 
;George Allen, March 2015
;
;Connects centerlines in the RivWidth v.0.4 output
; 
;Stand alone IDL program to be run after RivWidth v.0.4 runs.
;
;Inputs:
;RivWidth width image (and header file)
;RivWidth river mask image (and header file)
;
;Given an root directory, this program will run on all width and river mask files (including in subfolders) automatically.
;
;Outputs:
;Modified RivWidth width image:
;gaps filled (gapfill pixel DN=1)
;short spurs removed (e.g. river segments less than 100 px)
;
;Output files are written to a specified directory with a file name that contains the original w_image name.
;
;Routine:
;read input files
;calculate widths of river mask
;label regions of river mask
;fill in centerline gaps
; find centerline endpoints of width image
; subset around each centerline with a window size equal to the max width of the river mask
; if there is only one width centerline region, skip to next endpoint
; for each width centerline region other than the given endpoint region, find the closest pixel
; take a small subset (21 px) around each closest pixel and calculate the average width of the river mask centerline
; determine which is the best point for the endpoint to be connected to based on the width of each connecting points and their distance to the endpoint
; if there are no connecting points that have a high width to distance ratio, skip to the next endpoint
; subset around endpoint with a window size based on the distance between the endpoint and the connecting point
; label regions of the subset river mask. If the endpoint and connecting point exist in different regions then skip to the next endpoint
; connect endpoint to connection point with a single pixel line where DN=1
;remove duplicated lines by inverting width image and labeling regions
; fill regions less than 3 pixels (they will be thinned down later)
; for regions greater than 3 pixels, subset the islands
;   remove triple points from width centerline and label regions
;   if there are centerline pixels outside of the river mask, remove the region with the most centerline pixels outside the river mask
;   if there are no pixels outside the river mask (or multiple regions with an equal number of pixels outside the river mask), then remove the largest centerline region
; remove all triple point pixels that do not necessarily preserve connectivity
;remove any pixels that do not necessarily preserve connectivity
;remove any short centerline spurs that contain an endpoint
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; specify path locations:
rootDir = 'E:\GRWD\globalRun\rivWidth'
outDir = 'E:\GRWD\globalRun\shapefiles\gapless'

lastChar = strmid(rootDir, 0, 1, /reverse_offset)
if lastChar ne '\' or lastChar ne '/' then rootDir = rootDir + '\'
lastChar = strmid(outDir, 0, 1, /reverse_offset)
if lastChar ne '\' or lastChar ne '/' then outDir = outDir + '\'

; small subset window size (pixels):
wWinSize = 10
maxDim = 10

; window buffer size:
Bufr = 20

; define kernels:
kernel0 = [[1b,1,1],[1,0,1],[1,1,1]]
kernel1 = [[1b,1,1],[1,1,1],[1,1,1]]
kernel2 = [[2b,2,2],[2,2,2],[2,2,2]]
kernel3 = [[0b,2,0],[2,2,2],[0,2,0]]
kernel4 = [[0b,1,0],[1,0,1],[0,1,0]]
xInd = [-1,0,1,-1,0,1,-1,0,1]
yInd = [-1,-1,-1,0,0,0,1,1,1]

; get image paths:
rivPaths = file_search(rootDir, '*riv_class.img', /fold_case)
rivHdrsPaths = file_search(rootDir, '*riv_class.hdr', /fold_case)

wPaths = file_search(rootDir, '*w_image.img', /fold_case)
wHdrPaths = file_search(rootDir, '*w_image.hdr', /fold_case)

; get image name (might need to be modified depending on file naming convention):
fNames = strmid(wPaths, strlen(wPaths[0])-16, 4)

; check that the number of each file type is equal:
nWpaths = n_elements(wPaths)
if n_elements(rivPaths) ne nWpaths or n_elements(rivHdrsPaths) ne nWpaths $
  or n_elements(wHdrPaths) ne nWpaths then message, "Number of each file type is not equal!"

tm = systime(1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; run the gap filler routine on each mosaic image:  


asdf = [360]
for hh = 0, n_elements(asdf)-1 do begin ; 13, 154, 389, 624, 625, 638 ; nWpaths-1 do begin ;
TIC
h = asdf[hh]
print, ''
print, "h = " + strtrim(strtrim(string(h),1), 1) + "  ##########################  " + fNames(h) + string(tm)
print, "Opening file: " + wPaths[h] + "..."

envi_open_file, rivPaths[h], /no_interactive_query
envi_open_file, wPaths[h], /no_interactive_query

; read in river mask:
rivHdr = read_envi_hdr(rivHdrsPaths[h])
rImg = read_binary(rivPaths[h], data_dims=[rivHdr.samples, rivHdr.lines], data_type=1)

; read in width image:
wHdr = read_envi_hdr(wHdrPaths[h])
wImg = read_binary(wPaths[h], data_dims=[rivHdr.samples, rivHdr.lines], data_type=2)

; if there is no width centerlines, skip tile:
print,  n_elements(histogram(wImg))
if n_elements(histogram(wImg)) lt 2 then print, "No RivWidth centerline! Skipping tile..."
if n_elements(histogram(wImg)) lt 2 then continue

; find image dimensions:
width = wHdr.samples
height = wHdr.lines
offset = wHdr.header_offset

; get map info:
mapInfo = strsplit(wHdr.map_info, ",", /extract )
tp = mapInfo[3:4]

; find image resoluton and mean pixel spacing:
xRes = long(round(float(mapInfo[5])))
yRes = long(round(float(mapInfo[6])))
pxRes = (xRes + (2*xRes^2)^.5)/2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; calculate widths of river mask using thin and morph distance functions:
;if file_test(outDir+strtrim(string(h),1)+'rDist.tiff') eq 0 then begin
;  print, "Generating rDist image..."
  rDist = fix(morph_distance(rImg, neighbor_sampling=3))
  ;rThin = thin(rImg)
  ;rDist = rDist * rThin
;  write_tiff, outDir+strtrim(string(h),1)+'rDist.tiff', rDist, /short
;endif
;if file_test(outDir+strtrim(string(h),1)+'rDist.tiff') then $
;  rDist = read_image(outDir+strtrim(string(h),1)+'rDist.tiff')

; initial subset window size dependent on the max rDist width value contained in scene:
subWinSize = max(rDist)
;subWinSize = round(max(wImg)/pxRes) + 20

print, "subWindow size: " + strtrim(string(subWinSize),1) + " px"


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find and subset around width image centerline endpoints:
print, "Finding centerline endpoints..."

; label regions of wImg:
wReg = label_region(wImg, /all_neighbors, /uLong)

; find wImg endpoints (single neighbor pixels not on img edge):
wByte = bytarr(width, height)
wByte(where(wImg ne 0)) = 1
EPall = convol(wByte, kernel0, /edge_zero)
wByte(where(wByte eq 1 and EPall eq 1 and wReg ne 0)) = 2
EPall = array_indices(wByte, where(wByte eq 2))

; get endPt region values:
EPregVal = wReg[EPall[0, *], EPall[1, *]]
print, "Found " + strtrim(string(n_elements(EPregVal)),1) +" potential centerline gaps"

; label regions of wImg:
wReg = label_region(wImg, /all_neighbors, /uLong)

; subset images around end points:
subCoord = arraySubsetter(EPall, subWinSize, width, height)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; display endpoints (optional):
;imgSize =9
;
;x = wByte-wByte
;x[EPall[0,*], EPall[1,*]] = 255
;x = convol(x, kernel1)
;x = convol(x, kernel1)
;window, 0, xsize = width/imgSize, ysize = height/imgSize
;tvscl, reverse(congrid(x,width/imgSize,height/imgSize), 2), channel=1, /device
;tvscl, reverse(congrid(wByte,width/imgSize,height/imgSize), 2), channel=2
;tvscl, reverse(congrid(rImg,width/imgSize,height/imgSize), 2), channel=3
;txt = string([0:n_elements(EPall)/2-1])
;xyouts, EPall[0,*]/imgSize, (height-EPall[1,*])/imgSize, txt, alignment=1, charsize=.5, /device


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Connect endpoints

; For each endpoint, subset with a large window, find the closest pixels in each centerline region
; that is not the region of the endpoint, take a small subset of the river mask derived centerline
; and calculate the mean width. Find which centerline pixel is most likely to be the connecting point
; by balancing the connecting point width and distance from the endpoint. 

print, "Filling river centerline gaps... "

; loop through each endpoint:
for i = 14s, n_elements(subCoord[0,*])-1 do begin 

  ; initial large subset around endPts to search for nearby centerlines:
  subWr = wReg[subCoord[0,i]:subCoord[2,i]-1, subCoord[1,i]:subCoord[3,i]-1]
  
  ; if no other region(s) exists in sub image, skip to next endpoint:
  if n_elements(subWr(uniq(subWr, sort(subWr)))) lt 3 then continue
  
  ; For each region other than the end point region, find the closest pixel to end point:
  otherWind = where(subWr ne EPregVal(i) and subWr ne 0)
  ind = [0:n_elements(otherWind)-1]
  otherWreg = subWr(otherWind)
  
  ; sort nearby points by distance:
  otherWregXY = array_indices(subWr, otherWind)
  EPdistSq = (subCoord[4, i]-otherWregXY[0, *])^2 + (subCoord[5, i]-otherWregXY[1, *])^2
  distSort = sort(EPdistSq)
  sortReg = otherWreg(distSort)
  uniqReg = otherWreg(uniq(otherWreg, sort(otherWreg)))
  
  ; this could probably be done without the for loop using value_locate:
  CPraw = replicate(0, n_elements(uniqReg)) 
  for j = 0, n_elements(uniqReg)-1 do begin
    uniqRegLocs = where(uniqReg(j) eq sortReg) 
    CPraw(j) = uniqRegLocs(0)
  endfor
  
  CPind = otherWind(distSort(CPraw))
  CPdist = EPdistSq(distSort(CPraw))^.5
  subCP = array_indices(subWr, CPind)
  
  ; Display:
  subWr = wReg[subCoord[0,i]:subCoord[2,i]-1, subCoord[1,i]:subCoord[3,i]-1]
  ;subRdist = rDist[subCoord[0,i]:subCoord[2,i]-1, subCoord[1,i]:subCoord[3,i]-1]
  ;subWr[subCoord[4,i],subCoord[5,i]] = 9000
  ;subWr[subCP[0,*],subCP[1,*]] = 5000
  ;x = subWr & fName = 'subWr' & write_tiff, outDir+strtrim(string(h),1)+fName+'.tiff', x, /short  & envi_open_file, outDir+strtrim(string(h),1)+fName+'.tiff', /no_interactive_query
  ;x = wImg & fName = 'wImg' & write_tiff, outDir+strtrim(string(h),1)+fName+'.tiff', x, /short  & envi_open_file, outDir+strtrim(string(h),1)+fName+'.tiff', /no_interactive_query
  
  ; convert CPs coords from subImg to ImgIn and subset:
  CPs = [subCP[0,*]+subCoord[0,i], subCP[1,*]+subCoord[1,i]]
  sCoord = arraySubsetter(CPs, wWinSize, width, height)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; for each nearby centerline, find the max width of an area around the closest pixels:
  wDist = replicate(0, n_elements(CPind))
  for j = 0, n_elements(wDist)-1 do wDist(j) = max(rDist[sCoord[0,j]:sCoord[2,j]-1, sCoord[1,j]:sCoord[3,j]-1])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; for each nearby centerline, find the mean width of an area around the closest pixels: 
;  wDist = replicate(0, n_elements(CPind))
;  for j = 0, n_elements(wDist)-1 do begin
;    subRdist = rDist[sCoord[0,j]:sCoord[2,j]-1, sCoord[1,j]:sCoord[3,j]-1]
;    wDist(j) = mean(subRdist(where(subRdist gt 0)))
;  endfor
;  
;  ; if there are no nearby rivMaks derived centerlines, try using the rivWidth centerline data:
;  noCL = where(wDist eq 0, /null)
;  if noCL ne !null then begin
;    for j = 0, n_elements(noCL)-1 do begin
;      subWin = wImg[sCoord[0,noCL(j)]:sCoord[2,noCL(j)]-1, sCoord[1,noCL(j)]:sCoord[3,noCL(j)]-1]
;      wDist(noCL(j)) = mean(subWin(where(subWin gt 0)))/pxRes
;    endfor
;  endif
  
  ; find which of the nearby centerlines are closest relative to their river width:
  normDist = CPdist*1.3 - wDist - 100
  minDist = min(normDist)
  
  ; if wDist + constant > dist then nearby width, then skip to next endpoint:
  if minDist gt 0 then continue

  closestP = where(normDist eq minDist, /null)
  CP = CPs[*, closestP]
  
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; final subset based on the distance between CP and EP: 
  winSize = 2*round(CPdist(closestP))+100
  winSize = winSize(0)
  
  sCoord = arraySubsetter(EPall[*, i], winSize, width, height)
  
  ; get new XY coords of endPts and closestPts in the subset images:
  EP = sCoord[4:5]
  CP[0, where(sCoord[0] gt 0, /null)] = CP[0, 0] - EPall[0, i] + winSize
  CP[1, where(sCoord[1] gt 0, /null)] = CP[1, 0] - EPall[1, i] + winSize
  
  ; create final subset arrays:
  w = wImg[sCoord[0]:sCoord[2]-1, sCoord[1]:sCoord[3]-1]
  chan = rImg[sCoord[0]:sCoord[2]-1, sCoord[1]:sCoord[3]-1]
  
  ; determine if EP and CP locally span different river mask regions:
  rReg = label_region(chan, /all_neighbors)
  if rReg[EP[0, 0], EP[1, 0]] ne rReg[CP[0, 0], CP[1, 0]] then continue
  
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; connect centerlines with a 1 px line of DN = 1:
  
  ; draw a line between EP and CP:
  m = float(EP[1, 0] - CP[1, 0])/(EP[0, 0] - CP[0, 0])
  if m lt -1e5 or m gt 1e5 then m = 1e5
  b = EP[1, 0] - m*EP[0, 0]
  
  ; find whether the line is longer in the x or y dimension:
  if n_elements([EP[0, 0]:CP[0, 0]]) gt n_elements([EP[1, 0]:CP[1, 0]]) then begin
    xCon = [EP[0, 0]:CP[0, 0]]
    yCon = fix(round(m*xCon+b))
  endif else begin
    yCon = [EP[1, 0]:CP[1, 0]]
    xCon = fix(round((yCon-b)/m))
  endelse
  
  gapFill = bytarr(sCoord[6], sCoord[7])
  gapFill[xCon, yCon] = 1
  gapFill(where(w ne 0)) = 0
  
  ; add gapFill line back onto original width image:
  gapFill = gapFill + w
  wImg[sCoord[0]:sCoord[2]-1, sCoord[1]:sCoord[3]-1] = gapFill

endfor


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; remove duplicate fill lines by inverting the width image, labeling regions, and classifiying 
; the big and small islands. Fill

print, "Removing duplicated gap fill lines..."

; find islands between duplicated fill lines:
invReg = make_array(width, height, value=1, /byte)
invReg(where(wImg gt 0)) = 0
invReg = label_region(invReg, /uLong)
invRegHist = histogram(invReg)
uniqInvReg = [0s:n_elements(invRegHist)-1]

; if there are any centerline "islands," remove them:
if n_elements(invRegHist) ge 2 then begin
  ; set hist regions that touch the image boarder to zero: 
  edgeReg = [invReg[*, [1, height-2]], transpose(invReg[[1, width-2],*])]
  uniqEdgeReg = edgeReg(uniq(edgeReg, sort(edgeReg)))
  for i = 0, n_elements(uniqEdgeReg)-1 do invRegHist(where(uniqEdgeReg(i) eq uniqInvReg)) = 0
  
  ; categorize islands by their N pixels:
  smallReg = where(invRegHist gt 0 and invRegHist le 3, /null)
  bigReg = where(invRegHist gt 3, /null)
  
  ; fill in small regions (these will be thinned to a single pixel value):
  for i = 0, n_elements(smallReg)-1 do wImg(where(smallReg(i) eq invReg, /null)) = 1
  
  ; loop through each big island:
  for i = 0, n_elements(bigReg)-1 do begin
    ; determine coordinates of big islands:
    iRegXY = array_indices(invReg, where(bigReg(i) eq invReg))
    xMax = max(iRegXY[0, *], min=xMin)
    yMax = max(iRegXY[1, *], min=yMin)
  
    ; handle edges:
    xMin = xMin-bufr
    xMax = xMax+bufr
    yMin = yMin-bufr
    yMax = yMax+bufr
    xMin(where(xMin lt 0, /null)) = 0
    xMax(where(xMax gt width-1, /null)) = width
    yMin(where(yMin lt 0, /null)) = 0
    yMax(where(yMax gt height-1, /null)) = height
    dupWidth = xMax - xMin
    dupHeight = yMax - yMin
    
    ; subset to buffered region:
    rDup = rImg[xMin:xMax-1, yMin:yMax-1]
    wDup = wImg[xMin:xMax-1, yMin:yMax-1]

    ; find all triple points (three-neighbor pixels of wDup):
    wConv = bytarr(dupWidth, dupHeight)
    wConv(where(wDup gt 0)) = 2
    wConv = convol(wConv, kernel2, /edge_zero)
    tripPts = where(wConv ge 16 and wDup gt 0)
    
    ; remove three-neighbor pixels of wDup and label resulting regions:
    wDupRegs = wDup
    wDupRegs(tripPts) = 0
    wDupRegs = label_region(wDupRegs, /all_neighbors)
    
    ;  check if there are wDup pixels outside of riv mask:
    pxOutOfChan = where(wDupRegs ne 0 and rDup eq 0, /null)
    
    ; if there are pixels outside of riv mask then label their regions, 
    ; remove the region with the most pixels outside the region:
    if n_elements(pxOutOfChan) gt 0 then begin
      pxOutOfChanReg = wDupRegs(pxOutOfChan)
      pxOutOfChanRegHist = histogram(pxOutOfChanReg, min=0, max=max(pxOutOfChanReg))
      uniqPxOutOfChanReg = where(pxOutOfChanRegHist gt 0, /null)
      ; remove the region with the most pixels outside the rivMask:
      maxPxOutOfChanReg = where(pxOutOfChanRegHist eq max(pxOutOfChanRegHist))
      if n_elements(maxPxOutOfChanReg) eq 1 then wDup(where(wDupRegs eq maxPxOutOfChanReg(0))) = 0 $
      
      ; if there multiple regions with equal number of outside channel pixels,
      ; skip to the next step where the longer of the two regions is removed:
      else pxOutOfChan = !null
    endif
    
    ; if there are no pixels or multiple regions outside the river mask with an equal number of pixels, 
    ; remove the longest centerline of the loop from wDup:
    if n_elements(pxOutOfChan) eq 0 then begin
    
      uniqDupReg = wDupRegs(uniq(wDupRegs, sort(wDupRegs)))
      
      ; find which triple points that are connected to the main river centerlines:
      tripPtImg = bytarr(dupWidth, dupHeight)
      tripPtImg(tripPts) = 1
      tripPtReg = label_region(tripPtImg, /all_neighbors)
      tripPtDil = dilate(tripPtReg, kernel1, /gray)
      
      dupEdgeReg = [wDupRegs[*, [1, dupHeight-2]], transpose(wDupRegs[[1, dupWidth-2],*])]
      uniqDupEdgeReg = dupEdgeReg(uniq(dupEdgeReg, sort(dupEdgeReg)))
      dupEdgeReg = bytarr(dupWidth, dupHeight)
      for j = 1, n_elements(uniqDupEdgeReg)-1 do dupEdgeReg(where(uniqDupEdgeReg(j) eq wDupRegs)) = 1
      dupEdgeReg(where(tripPtDil gt 0 and dupEdgeReg gt 0)) = 2
      mainTripPts = tripPtDil(where(dupEdgeReg eq 2))
      
      ; delete the larger loop region from wDup:
      mainTripRm = wDup
      for j = 0, n_elements(mainTripPts)-1 do mainTripRm(where(mainTripPts(j) eq tripPtReg)) = 0
      mainTripRmReg = label_region(mainTripRm, /all_neighbors)
      mainTripRmReg(where(mainTripRmReg gt 0 and dupEdgeReg gt 0)) = 0
      mainTripRmRegHist = histogram(mainTripRmReg)
      mainTripRmRegHistInd = [1:n_elements(mainTripRmRegHist)-1]
      mainTripRmRegHist = mainTripRmRegHist[mainTripRmRegHistInd]
      sortedRegs = reverse(sort(mainTripRmRegHist))
      wDup(where(mainTripRmRegHistInd(sortedRegs(0)) eq mainTripRmReg)) = 0
    
    endif

    ; insert wImg back into wImg image: 
    wImg[xMin:xMax-1, yMin:yMax-1] = long(wDup[0:dupWidth-1, 0:dupHeight-1])
    
    ; find the maximum island dimension to be used to subset triple points later:
    maxDim = max([maxDim, dupWidth, dupHeight])
    
  endfor
  
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Centerline clean up:

; remove pixels that do not preserve connectivity or short spurs:

; identify triple points px and determine if they are needed to preserve centerline connectivity:
print, "Removing nonessentail triple pt pixels..."

replicate_inplace, wByte, 0
wByte(where(wImg gt 0)) = 2
wByte = convol(wByte, kernel2, /edge_zero)
tripPts = where(wByte ge 16 and wImg ne 0)
if tripPts[0] ne -1 then begin
  tripXY = array_indices(wByte, tripPts) 
  tsC = arraySubsetter(tripXY, maxDim+10, width, height)
  wVals = wImg(tripPts)
  for i = 0, n_elements(tripPts)-1 do begin
    sub = wImg[tsC[0,i]:tsC[2,i]-1, tsC[1,i]:tsC[3,i]-1]
    nReg = max(label_region(sub, /all_neighbors))
    sub[tsC[4,i], tsC[5,i]] = 0
    if max(label_region(sub, /all_neighbors)) ne nReg then sub[tsC[4,i], tsC[5,i]] = wVals(i)
    wImg[tsC[0,i]:tsC[2,i]-1, tsC[1,i]:tsC[3,i]-1] = sub
  endfor
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; remove any uncessary pixels:
; loop through each centerline px, subset with a small window, and remove if the px is not 
; necessarily needed for connectivity. Exclude endpoints from the list:
print, "Removing all nonessentail pixels..."
replicate_inplace, wByte, 0
wByte(where(wImg gt 0)) = 1
EPall = convol(wByte, kernel0, /edge_zero)
EPall = where(wByte eq 1 and EPall eq 1, /null)
wByte(EPall) = 2
cLinePx = where(wImg gt 0 and wByte ne 2, /null)
cLineXY = array_indices(wImg, cLinePx) 
sC = arraySubsetter(cLineXY, wWinSize, width, height)
wVals = wImg(cLinePx)
for i = 0, n_elements(cLinePx)-1 do begin
  sub = wImg[sC[0,i]:sC[2,i]-1, sC[1,i]:sC[3,i]-1]
  nReg = max(label_region(sub, /all_neighbors))
  sub[sC[4,i], sC[5,i]] = 0
  if max(label_region(sub, /all_neighbors)) ne nReg then sub[sC[4,i], sC[5,i]] = wVals(i)
  wImg[sC[0,i]:sC[2,i]-1, sC[1,i]:sC[3,i]-1] = sub
endfor 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; remove short (e.g. <10 km) centerline spurs by removing triple points, and counting pixels of 
;  regions that touch an endpoint:
print, "Removing short centerline spurs..."
replicate_inplace, wByte, 0
wByte(where(wImg gt 0)) = 2
conv = convol(wByte, kernel2, /edge_zero)
tripPts = where(conv ge 16 and wImg ne 0)
if tripPts[0] ne -1 then begin
  wByte(tripPts) = 0
  regs = label_region(wByte, /all_neighbors, /uLong)
  replicate_inplace, wByte, 0
  wByte(EPall) = 1
  wByte = dilate(wByte, kernel1)

  ; remove short spurs (small regions) that contain endpoints:
  minSpurSize = 333 ;n px
  EPregHist = histogram(regs(where(regs gt 0 and wByte eq 1)), min=0, max=max(regs))
  regHist = histogram(regs)
  histInd = [0:n_elements(regHist)-1]
  nPixEndPtReg = regHist*(EPregHist gt 0)
  spurReg = histInd(where(nPixEndPtReg lt minSpurSize and nPixEndPtReg gt 0, /null))
  for i = 0, n_elements(spurReg)-1 do wImg(where(regs eq spurReg(i))) = 0
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; remove any uncessary pixels again:
; loop through each centerline px, subset with a small window, and remove if the px is not
; necessarily needed for connectivity. Exclude endpoints from the list:
print, "Removing all nonessentail pixels..."
replicate_inplace, wByte, 0
wByte(where(wImg gt 0)) = 1
EPall = convol(wByte, kernel0, /edge_zero)
EPall = where(wByte eq 1 and EPall eq 1)
wByte(EPall) = 2
cLinePx = where(wImg gt 0 and wByte ne 2)
cLineXY = array_indices(wImg, cLinePx)
sC = arraySubsetter(cLineXY, wWinSize, width, height)
wVals = wImg(cLinePx)
for i = 0, n_elements(cLinePx)-1 do begin
  sub = wImg[sC[0,i]:sC[2,i]-1, sC[1,i]:sC[3,i]-1]
  nReg = max(label_region(sub, /all_neighbors))
  sub[sC[4,i], sC[5,i]] = 0
  if max(label_region(sub, /all_neighbors)) ne nReg then sub[sC[4,i], sC[5,i]] = wVals(i)
  wImg[sC[0,i]:sC[2,i]-1, sC[1,i]:sC[3,i]-1] = long(sub)
endfor


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; remove any very small regions:

regs = label_region(wImg, /all_neighbors, /uLong)
regHist = histogram(regs(where(regs gt 0)), min=0, max=max(regs))
smallReg = where(regHist lt 10 and regHist gt 0, /null)
if smallReg ne !null then print, "Removing very small regions..."
for i = 0, n_elements(smallReg)-1 do wImg(where(regs eq smallReg(i), /null)) = 0


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; write out connected width image:
;x = wImg & fName = 'gapFilled'
;x = wImg & fName = 'gapFilled' & write_tiff, outDir+strtrim(string(h),1)+fNames(h)+'_'+fName+'.tif', x, /long 
;envi_open_file, outDir+strtrim(string(h),1)+fNames(h)+'_'+fName+'.tif', /no_interactive_query

TOC


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; export as a polyline shapefile:

print, 'Organizing pixels into vectors..."

TIC

; create new shapefile object and define the entitiy type to line:
sfName = outDir+fNames(h)
print, sfName
outshape = obj_new('IDLffShape', sfName+'.shp',  /UPDATE, ENTITY_TYPE=3)

; set the attribute definitions for the new Shapefile
outshape->addAttribute, 'easting', 3, 12, precision=0
outshape->addAttribute, 'northing', 3, 12, precision=0
outshape->addAttribute, 'width', 3, 7, precision=0
outshape->addAttribute, 'nchannels', 3, 3, precision=0
outshape->addAttribute, 'segmentID', 3, 7, precision=0
outshape->addAttribute, 'segmentInd', 3, 7, precision=0

entcounter = 0l

; organize river segments and triple points:

; Once again find triple points: 
replicate_inplace, wByte, 0
wByte(where(wImg ne 0)) = 2

conv = convol(wByte, kernel2, /edge_zero)
tripPts = where(conv ge 16 and wImg ne 0, /null)

; set triple points to zero and label regions:
wByte(tripPts) = 0
regs = label_region(wByte, /all_neighbors, /uLong)

; sort river segments from large to small (probably not be necessary):
regHist = histogram(regs)
sortReg = reverse(sort(regHist))
sortReg = sortReg[1:n_elements(regHist)-1]
;if tripPts ne !null then begin
;  tripXY = array_indices(wImg, tripPts)
;  tripSC = arraySubsetter(tripXY, 1, width, height)
;endif


; label triple point regions:
replicate_inplace, wByte, 0
wByte(tripPts) = 1


; convert wImg from int to long int (for very wide widths):
; This is funky. One solution could be to devide widths by resolution of imagery 
; to avoid having record very big numbers. 
wImg = long(wImg)
wImg[where(wImg lt 0)] = wImg[where(wImg lt 0)]+65536

; speed things up by subseting each river region before determining order:
; fixme: consider spliting up long river regions to possibly speed up this process
bufr = 5

for i = 0, n_elements(sortReg)-1 do begin
  
  ; find XY of river region:
  segXY = array_indices(regs, where(sortReg(i) eq regs))
  
  ; subset around this shape with a buffer:
  xMax = max(segXY[0, *], min=xMin)
  yMax = max(segXY[1, *], min=yMin)

  ; handle edges:
  xMin = xMin-bufr
  xMax = xMax+bufr
  yMin = yMin-bufr
  yMax = yMax+bufr
  xMin(where(xMin lt 0, /null)) = 0
  xMax(where(xMax gt width-1, /null)) = width
  yMin(where(yMin lt 0, /null)) = 0
  yMax(where(yMax gt height-1, /null)) = height
  dupWidth = xMax - xMin
  dupHeight = yMax - yMin
  
  ; find segment subset XY:
  sXY = segXY
  sXY[0, *] = segXY[0, *] - xMin
  sXY[1, *] = segXY[1, *] - yMin
  
  ; subset labeled segment and triple triple points:
  sSeg = bytarr(dupWidth, dupHeight)
  sSeg[sXY[0, *], sXY[1, *]] = 1
  
  sTrip = wByte[xMin:xMax-1, yMin:yMax-1]
  tripReg = label_region(sTrip, /all_neighbors, /uLong)
  tripHist = histogram(tripReg, min=0, max=max(tripReg))
  
  ; identify the two segment endpoints:
  sEP = convol(sSeg, kernel0, /edge_zero)
  sEP = where(sSeg eq 1 and sEP eq 1, /null)
  
  ; in the case of a single pixel region, don't create a shapefile entity:
  if sEP eq !null then continue 
  
  ; check if there are any triple point pixels adjacent to the endpoints:
  sEPXY = array_indices(sSeg, sEP)
  sEPSC = arraySubsetter(sEPXY, 1, dupWidth, dupHeight)

  for j = 0, n_elements(sEP)-1 do begin 
    sEPsub = tripReg[sEPSC[0,j]:sEPSC[2,j], sEPSC[1,j]:sEPSC[3,j]]
    tRegMatchInd = where(sEPsub ne 0, /null)
    tRegMatch = sEPsub[tRegMatchInd]
    tripCheck = tripHist[tRegMatch]

    ; if the triple point px is a single pixel, add it to the segment and update endpoint list:
    if tripCheck eq 1 then begin
      tripMatchInd = where(tRegMatch(0) eq tripReg)
      sSeg[tripMatchInd] = 1
      sEP[j] = tripMatchInd
      
    endif
    
    ; if the triple point pixel:
    if tripCheck eq 3 then begin
      ; triple point XY closest to segment endpoint:
      sTripMatchXY = [sEPXY[0, j] + xInd(tRegMatchInd), sEPXY[1, j] + yInd(tRegMatchInd)]
      
      ; update endpoint list and add this pixel to the segment:
      sEP[j] = sTripMatchXY[1]*dupWidth + sTripMatchXY[0]
      sSeg[sEP[j]] = 1
      
      ; determine if the pixel is a double or triple endpoint location:
      sTSC = arraySubsetter(sTripMatchXY, 1, dupWidth, dupHeight)
      otherTripPx = tripReg[sTSC[0]:sTSC[2], sTSC[1]:sTSC[3]]
      otherTripPxInd = where(kernel4*otherTripPx ne 0, /null)
      
      ; if the pixel is a double pixel, update the endpoint location and add the triple pixel to the segment:
      if n_elements(otherTripPxInd) eq 1 then begin
        tTripXY = [sTripMatchXY[0] + xInd(otherTripPxInd), sTripMatchXY[1] + yInd(otherTripPxInd)]
        sEP[j] = tTripXY[1]*dupWidth + tTripXY[0]
        sSeg[sEP[j]] = 1
      endif
      
    endif
    
  endfor

  
  ; start at the first endpoint and grow region along line until second endpoint is reached:
  crawler = sSeg*2
  crawler(sEP(0)) = 4

  sXY = array_indices(sSeg, where(sSeg ne 0))
  cL = sXY
  
  for j = 0, n_elements(cL)/2-1 do begin
    cL[*, j] = array_indices(crawler, where(crawler eq 4))
    crawler[cL[0,j]-1:cL[0,j]+1, cL[1,j]-1:cL[1,j]+1] = crawler[cL[0,j]-1:cL[0,j]+1, cL[1,j]-1:cL[1,j]+1]*2   
  endfor
  
  ; add vertices and attribute data to shapefile:
  x = cL[0, *] + xMin
  y = cL[1, *] + yMin
  utmX = TP[0] + x * xRes + xRes/2
  utmY = TP[1] - y * yRes - yRes/2
  
  ; beacuse each segment bridges two pixels, take the mean width of the two pixels:
  N = n_elements(x)
  segW = (wImg[x[0:N-2], y[0:N-2]] + wImg[x[1:N-1], y[1:N-1]])/2
  
  ; any segment that touches a connection pixel (px=1), set segment value equal to 1: 
  cInd = where(wImg[x[0:N-1], y[0:N-1]] eq 1, /null)
  if cInd ne !null then begin
    cInd = [cInd, cInd-1]
    cInd = cInd(uniq(cInd, sort(cInd)))
    segW(cInd(where(cInd ge 0 and cInd le N-1, /null))) = 1
  endif
  
  ; one line segments between two pixels:
  for j = 0, N-2 do begin

    ; create structure for new entity:
    entSeg = {IDL_SHAPE_ENTITY}
    entSeg.Shape_type = 3
    entSeg.Bounds = [min(utmX[j:j+1]), min(utmY[j:j+1]), 0, 0, max(utmX[j:j+1]), max(utmY[j:j+1]), 0, 0]
    entSeg.n_vertices = 2
    entSeg.vertices = ptr_new([utmX[j], utmY[j], utmX[j+1], utmY[j+1]])
  
    ; add data to attribute table:
    attrSeg = outshape->getAttributes(/attribute_structure)
    attrSeg.attribute_0 =  utmX[j]
    attrSeg.attribute_1 = utmY[j]
    attrSeg.attribute_2 = segW[j]
    ; fixme need to make a channel image to do this. If properly integrated with RivWidth, this will not be necessary: 
    attrSeg.attribute_3 = 0
    attrSeg.attribute_4 = i
    attrSeg.attribute_5 = j
    
    outshape->putEntity, entSeg
    outshape->setAttributes, entcounter, attrSeg
    entcounter++

  endfor
  
endfor


; create prj file for shapefile:
shpPrjOut = sfName + '.prj'
openw, lun, shpPrjOut, /get_lun
printf, lun, wHdr.coord_system
close, lun
free_lun, lun

print, "Completed mosaic: " + fNames(h)

obj_destroy, outshape

TOC

print, "Finished in (seconds): ", systime(1) - tm

endfor

end





















;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; sub functions:

function read_envi_hdr, filename
  ; define header structure
  header = { description:"", $
    samples:-1l, $
    lines:-1l, $
    bands:-1l, $
    header_offset:-1l, $
    file_type:"", $
    data_type:-1l, $
    interleave:"", $
    sensor_type:"", $
    byte_order:-1l, $
    wavelength_units:"", $
    z_plot_range:[-1d,-1d], $
    z_plot_titles:["",""], $
    map_info:"", $
    coord_system:"", $
    band_names:ptr_new(), $
    wavelength:ptr_new() $
  }

  openr, lun, filename, /get_lun

  ; read the first record and determine if this is a valid
  ; envi header file
  str = ""
  readf, lun, str
  str = strtrim( strcompress( str ), 2 )
  if ( str ne "ENVI" ) then begin
    message, "envi header file has an invalid format", /continue
    return, -1
  endif

  ; read each subsequent record until the end of file is reached
  while not eof( lun ) do begin

    readf, lun, str

    ; if the current record does not have zero length, proceed to
    ; parse the information for this record
    if ( strlen( str ) gt 0 ) then begin

      ; determine if the current record is the beginning of a new
      ; name/value pair or if it is the continuation of a multiple
      ; line value and parse the pair appropriately
      str = strcompress( str )
      equalposition = strpos( str, "=" )
      if ( equalposition gt 0 ) then begin
        name = strtrim( strmid( str, 0, equalposition ), 2 )
        value = strtrim( strmid( str, equalposition+1 ), 2 )
        multilinevalue = strtrim( value, 2 )
      endif else begin
        multilinevalue = multilinevalue + " " + strtrim( str, 2 )
      endelse

      ; if the value is defined across multiple lines, see if both the
      ; beginning and ending delimiters are present: if not, concatenate
      ; the current line with the previous line and continue reading data,
      ; otherwise, strip the delimiters and proceed to parse the completed
      ; value
      if ( strpos( multilinevalue, "{" ) ne -1 ) and ( strpos( multilinevalue, "}" ) ne -1 ) then begin
        multilinecomplete = 1
        multilinevalue = strtrim( multilinevalue, 2 )
        multilinevalue = strmid( multilinevalue, 1, strlen( multilinevalue )-2 )
      endif else begin
        multilinecomplete = 0
      endelse

      ; parse the name/value pair
      case strlowcase( name ) of
        "description": if multilinecomplete then header.description = multilinevalue
        "samples": header.samples = long( value )
        "lines": header.lines = long( value )
        "bands": header.bands = long( value )
        "header offset": header.header_offset = long( value )
        "file type": header.file_type = value
        "data type": header.data_type = long( value )
        "interleave": header.interleave = strlowcase( value )
        "sensor type": header.sensor_type = value
        "byte order": header.byte_order = long( value )
        "wavelength units": header.wavelength_units = value
        "z plot range": if multilinecomplete then header.z_plot_range = double( strsplit( multilinevalue, ",", /extract ) )
        "z plot titles": if multilinecomplete then header.z_plot_titles = strsplit( multilinevalue, ",", /extract )
        "band names": if multilinecomplete then header.band_names = ptr_new( strsplit( multilinevalue, ",", /extract ) )
        "wavelength": if multilinecomplete then header.wavelength = ptr_new( double( strsplit( multilinevalue, ",", /extract ) ) )
        "map info":  if multilinecomplete then header.map_info = multilinevalue
        "coordinate system string": if multilinecomplete then header.coord_system = multilinevalue

        else:
      endcase

    endif

  endwhile

  ; release the current logical unit number
  free_lun, lun

  ; return the filled header structure to the calling routine
  return, header

end


function arraySubsetter, XY, wWin, imgWidth, imgHeight
  ; clips arrays into a square subet including at edges of arrays. 
  ; input of a center point to subset around, width of array (width/2), big array width and height
  ; ouputs coords of subset and new coords of central point (in the subset coordinates), the subset width and height
  p = XY
  xMin = p[0, *] - wWin
  yMin = p[1, *] - wWin
  xMax = p[0, *] + wWin
  yMax = p[1, *] + wWin
  
  p[0, where(xMin gt 0, /null)] = wWin
  p[1, where(yMin gt 0, /null)] = wWin
  
  xMin(where(xMin lt 0, /null)) = 0
  yMin(where(yMin lt 0, /null)) = 0
  xMax(where(xMax gt imgWidth, /null)) = imgWidth
  yMax(where(yMax gt imgHeight, /null)) = imgHeight
  
  subWidth = xMax - xMin
  subHeight = yMax - yMin
  
  subCoords = [xMin, yMin, xMax, yMax, p[0,*], p[1,*], subWidth, subHeight]
  return, subCoords
end






