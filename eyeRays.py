# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.lines import Line2D
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.transforms as transforms

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = "8"

#######################################################################################
def myWedge(ax, cx=0.5, cy=0.5, rad=0.3, theta1=60, theta2=120,
                ex=0.1, ey=0.5, **kwargs):
  pathArc = Path.arc(theta1, theta2) # arc of unit circle
  transRad = transforms.Affine2D().scale(rad)
  transCenter = transforms.Affine2D().translate(cx, cy)
  transArc = transRad + transCenter
  pathArc = pathArc.transformed(transArc)
  arc_start = pathArc.vertices[0]
  arc_end = pathArc.vertices[-1]
  vertices = pathArc.vertices
  codes = pathArc.codes
  vertices = np.append(vertices, [[ex, ey], arc_start], 0)
  codes = np.append(codes, [Path.LINETO]*2)
  pathpatch = PathPatch(Path(vertices, codes), transform=ax.transData, **kwargs)
  ax.add_patch(pathpatch)
################################
def drawEyeRays(ax, params):
  ################
  def cosine_theorem(a, b, c):
    cosA = (a*a + b*b - c*c)/(2*a*b)
    return np.arccos(cosA)
  ################

  ax.set_xlim(-15, 15)
  ax.set_ylim(-12.8, 22.5)
  ax.set_aspect('equal', 'box')
  ax.set_xticks([])
  ax.set_yticks([])
  #ax.axis('off')
  
  Deye = params['Deye']
  Reye = Deye/2.0
  Rcor = params['Rcor']
  ACD = params['ACD']
  Ncor = params['Ncor']

  OccluderDiam = params['OccluderDiam']
  PupilDiam = params['PupilDiam']

  EccDeg = params['EccDeg']
  Ecc    = np.deg2rad(EccDeg) # Eccentricity (rad)
  ThetaRad = 0.5*np.pi - Ecc

  DeltaRays = params['DeltaRays']
  MaxNumRays = params['MaxNumRays']
  MainRayLength = params['MainRayLength']

  a = Rcor - ACD
  # Note that:
  # b = Rcor*np.sin(phi2)
  # b = Reye*np.sin(phi1)
  phiCor = np.arccos(a/Rcor)
  phiEye = np.arcsin(Rcor/Reye*np.sin(phiCor))
  #print('phiEye='+str(phiEye)+ '('+ str(np.rad2deg(phiEye)) + ')')
  #print('phiCor='+str(phiCor)+ '('+ str(np.rad2deg(phiCor)) + ')')

  deltaReyeRcor = Reye*np.cos(phiEye) - a

  ### Eye Ball
  angle=90
  eyeBall = patches.Arc(xy=(0, 0), width=2*Reye, height=2*Reye, angle=angle,
                           theta1=+np.rad2deg(phiEye), theta2=-np.rad2deg(phiEye))
  ax.add_patch(eyeBall)
  ### Cornea
  cornea   = patches.Arc(xy=(0, deltaReyeRcor), width=2*Rcor, height=2*Rcor, angle=angle,
                           theta1=-np.rad2deg(phiCor), theta2=+np.rad2deg(phiCor))
  ax.add_patch(cornea)

  ### Iris
  yIris = deltaReyeRcor + Rcor*np.cos(phiCor)
  xIris = Rcor*np.sin(phiCor)
  xPupil = PupilDiam/2.0
  ax.plot([-xIris, -xPupil], [yIris, yIris], lw=2, color='black')
  ax.plot([xPupil, xIris], [yIris, yIris], lw=2, color='black')

  #Occluder
  deltaRocc = 0.4
  Rocc = Rcor + deltaRocc
  phiOccluder = 0.5*OccluderDiam/Rocc
  #print('phiOccluder_deg='+str(np.rad2deg(phiOccluder)))
  occluder = patches.Arc(xy=(0, deltaReyeRcor), width=2*Rocc, height=2*Rocc, angle=angle,
                           theta1=-np.rad2deg(phiOccluder), theta2=+np.rad2deg(phiOccluder), lw=1, zorder=10)
  ax.add_patch(occluder)

  #########################
  #Rays

  rays = []
  NumRays = MaxNumRays
  for i in range(0, MaxNumRays+1):
    ## segment0 - from infinity to cornea
    #mainRay passes through the right edge of cornea
    if (i==0):
      if (Ecc < 0):
        phiCor = -phiCor
        DeltaRays = -DeltaRays
      endX0 = Rcor*np.sin(phiCor)
      endY0 = Rcor*np.cos(phiCor) + deltaReyeRcor
      startX0 = endX0 + MainRayLength*np.sin(Ecc)
      startY0 = endY0 + MainRayLength*np.cos(Ecc)

    startX = startX0 - i*DeltaRays*np.cos(Ecc) # note cos!
    startY = startY0 + i*DeltaRays*np.sin(Ecc) # note sin!
    # helper ray passing through the center of the cornea circle
    # D = D1 + D2
    # D2 = Rocc*sin(phiIncident)
    D = i*DeltaRays
    D1 = Rcor*np.sin(phiCor - Ecc)
    #D2 = {Rcor|Rocc}*sinPhiIncident    
    D2 = D - D1
    sinPhiIncidentOcc = D2/Rocc
    sinPhiIncidentCor = D2/Rcor
    if (abs(sinPhiIncidentOcc)>=1.0 or startY < yIris):
      NumRays = i
      break
    if (abs(sinPhiIncidentCor)>=1.0 or startY < yIris):
      NumRays = i
      break
    phiIncidentOcc = np.arcsin(sinPhiIncidentOcc)
    phiIncidentCor = np.arcsin(sinPhiIncidentCor)
    phiRayOcc = Ecc - phiIncidentOcc # incident angle
    phiRayCor = Ecc - phiIncidentCor # incident angle
    ## passing Occluder or not?
    rayStoppedByOccluder = False
    if (phiRayOcc >= phiOccluder or phiRayOcc <= -phiOccluder):
      drawArrow=False
    else:
      rayStoppedByOccluder = True
      drawArrow=True
    #print('i='+ str(i)+'; phiRayOcc='+ str(np.rad2deg(phiRayOcc))+'; phiRayCor='+ str(np.rad2deg(phiRayCor)))
    if(i==0):
      endX = endX0
      endY = endY0
    else:
      if (rayStoppedByOccluder):
        endX = Rocc*np.sin(phiRayOcc)
        endY = Rocc*np.cos(phiRayOcc) + deltaReyeRcor
      else:
        endX = Rcor*np.sin(phiRayCor)
        endY = Rcor*np.cos(phiRayCor) + deltaReyeRcor
    if (endY < yIris):
      NumRays = i
      break
     
    segment0 = {'startX': startX, 'startY': startY,
                'endX': endX, 'endY': endY,
                'drawArrow':drawArrow, 'rayEcc':Ecc, 'drawRay':True}
    ##################################################
    ### segment1: cornea -> iris or back side of eye
    phiRefr = np.arcsin(np.sin(phiIncidentCor)/Ncor)
    #print('i='+ str(i)+'; phiIncident='+str(np.rad2deg(phiIncident))+'; phiRefr='+str(np.rad2deg(phiRefr)))
    eccRefr = phiRayCor + phiRefr
    #print('i='+str(i)+'; eccRefr='+str(np.rad2deg(eccRefr)))
    startX = endX
    startY = endY
    drawRay = True
    if (i==0):
      drawRay = False
    if (i == NumRays):
      drawRay = False
    ########### does ray stopped by iris?
    rayStoppedByIris = False
    L = (startY - yIris)/np.cos(eccRefr)
    x = startX - L*np.sin(eccRefr)
    if(x > xPupil or x < -xPupil):
      rayStoppedByIris = True
      endX = x       
      endY = yIris
      drawArrow= True
    else:
      d1 = Rcor*np.sin(phiRefr)
      d2 = deltaReyeRcor * np.sin(eccRefr)
      D = d1 + d2
      phiRayEye = np.arcsin(D/Reye) + eccRefr # relative to Reye
      #print('i='+str(i)+'; phiRayEye='+str(np.rad2deg(phiRayEye)))
      endX = -Reye*np.sin(phiRayEye)
      endY = -Reye*np.cos(phiRayEye)
      drawArrow = False
    segment1 = {'startX': startX, 'startY': startY,
                'endX': endX, 'endY': endY,
                'drawArrow':drawArrow, 'rayEcc':eccRefr,
                'drawRay':drawRay,
                'rayStoppedByOccluder':rayStoppedByOccluder,
                'rayStoppedByIris':rayStoppedByIris}
    segments = [segment0, segment1]
    rays.append(segments)
    #print('ray['+str(i)+']='+str(rays[i]))
  
  ####################################################
  ######## plot all

  #plot cornea center
  #ax1.plot(0, deltaReyeRcor, marker='o', markersize=1, color='black')

  #plot eye center
  ax.plot(0, 0, marker='o', markersize=3, color='black')

  #plot vertical dashed line
  ax.plot([0, 0], [-13, +22], color='black', ls='--', lw=0.5)

  #set drawRay = False for intermediate 'gray' rays
  #only first and last gray ray should be drawn
  grayArea = False
  prevGrayArea = False
  segment = 1
  firstGrayRay = -1
  lastGrayRay = -1
  for i in range(0, NumRays):
    drawRay=rays[i][segment]['drawRay']
    rayStoppedByOccluder = rays[i][segment]['rayStoppedByOccluder']
    rayStoppedByIris = rays[i][segment]['rayStoppedByIris']
    if (rayStoppedByOccluder and not rayStoppedByIris):
      if (not grayArea):
        grayArea = True
        drawRay = True
        firstGrayRay = i
      else:
        drawRay = False
    if (rayStoppedByOccluder and rayStoppedByIris):
      if (prevGrayArea):
        drawRay = True
        rays[i-1][segment]['drawRay'] = drawRay
        lastGrayRay = i-1
        grayArea = False
    prevGrayArea = grayArea
    rays[i][segment]['drawRay'] = drawRay

  #print('firstGrayRay='+str(firstGrayRay)+'; lastGrayRay='+str(lastGrayRay))
  #######
  rayWidth = 0.001
  for i in range(0, NumRays):
    for segment in range(0, 2):
      #print('ray['+str(i)+']['+str(segment)+']='+str(rays[i][segment]))
      startX = rays[i][segment]['startX']
      startY = rays[i][segment]['startY']
      endX = rays[i][segment]['endX']
      endY = rays[i][segment]['endY']
      rayEcc = rays[i][segment]['rayEcc']
      drawArrow=rays[i][segment]['drawArrow']
      drawRay=rays[i][segment]['drawRay']
      linestyle = '-'
      if (segment == 1):
        rayStoppedByOccluder = rays[i][segment]['rayStoppedByOccluder']
        rayStoppedByIris = rays[i][segment]['rayStoppedByIris']
        if (rayStoppedByOccluder and not rayStoppedByIris):
          #linestyle = 'dotted'
          linestyle=(0, (5,10))
          drawArrow = False
        if (rayStoppedByOccluder and i > lastGrayRay):
          drawRay = False
      if (drawArrow):
  #     head_width=3*rayWidth
        head_width=0.1
        rayLength = np.sqrt((endX - startX)**2 + (endY - startY)**2) - 0.4
        endX = startX - rayLength*np.sin(rayEcc)
        endY = startY - rayLength*np.cos(rayEcc)
      else:
        head_width=0.0
      #print('ray['+str(i)+']['+str(segment)+']='+str(rays[i][segment])+ ' ls='+str(linestyle))
      if (drawRay):
        dx = endX - startX
        dy = endY - startY
        ax.arrow(startX, startY, dx, dy,
                length_includes_head=True, width=rayWidth, head_width=head_width, lw=0.5, ls=linestyle, color='blue')
  ## draw gray area
  if (firstGrayRay >=0 and lastGrayRay >=0):
    segment = 1
    startX1 = rays[firstGrayRay][segment]['startX']
    startY1 = rays[firstGrayRay][segment]['startY']
    endX1 = rays[firstGrayRay][segment]['endX']
    endY1 = rays[firstGrayRay][segment]['endY']
    rayEcc1 = rays[firstGrayRay][segment]['rayEcc']

    theta1 = np.rad2deg(np.arccos(startX1/Rcor))
  
    startX2 = rays[lastGrayRay][segment]['startX']
    startY2 = rays[lastGrayRay][segment]['startY']
    rayEcc2 = rays[lastGrayRay][segment]['rayEcc']
    theta2 = np.rad2deg(np.arccos(startX2/Rcor))
  
    rayLength1 = np.sqrt((endX1 - startX1)**2 + (endY1 - startY1)**2) - 0.2

#    myWedge(ax, cx=0, cy=deltaReyeRcor, rad=Rcor, theta1=theta1, theta2=theta2, ex=endX1, ey=endY1, color='gray')
    myWedge(ax, cx=0, cy=deltaReyeRcor, rad=Rocc, theta1=theta1, theta2=theta2, ex=endX1, ey=endY1, color='gray')
##########################################################################  
### END drawEyeRays ######################################################

if __name__ == "__main__":
  fig = plt.figure(figsize=(4.2, 4.2))
  fig.patch.set_visible(False)

  rect1  = [0.05, 0.05, 0.95, 0.95] #[left, bottom, width, height]
  ax1 = fig.add_axes(rect1)

  #PARAMETERS
  params = {}

  params['Deye'] = 24.5 #Eye diameter (mm)
  params['Rcor'] = 7.8  #Corneal curvature radius (mm)
  params['ACD']  = 3.1  #Anterior Chamber Depth (mm)
  params['Ncor'] = 1.38

  params['OccluderDiam'] = 5
  params['PupilDiam'] = 3

  params['EccDeg'] = 50 # Eccentricity (degrees)
  EccDeg = params['EccDeg']
  Ecc    = np.deg2rad(EccDeg) # Eccentricity (rad)
  ThetaRad = 0.5*np.pi - Ecc

  params['DeltaRays'] = 0.52
  params['MaxNumRays'] = 30
  params['MainRayLength'] = 7

  drawEyeRays(ax1, params)
  ax1.text(-15, 20, 'Eccentricity '+str(params['EccDeg'])+'°')
  ax1.text(-15, 18, 'Occluder '+str(params['OccluderDiam'])+' mm')
  ax1.text(-15, 16, 'Pupil '+str(params['PupilDiam'])+' mm')

#  plt.savefig('eyeRays_Ecc'+str(params['EccDeg'])+'_Occ'+str(params['OccluderDiam'])+'_Pupil'+str(params['PupilDiam'])+'.png', dpi=600)
  plt.show()

