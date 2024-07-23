
if True:

  umax = 40 # u surface divisions
  vmax = 20 # v surface divisions
  npts = (umax+1)*(vmax+1)

  #surface1
  s1x = [[0 for v in range(vmax+1)] for u in range(umax+1)] 
  s1y = [[0 for v in range(vmax+1)] for u in range(umax+1)] 
  s1z = [[0 for v in range(vmax+1)] for u in range(umax+1)] 
  #surface2
  s2x = [[0 for v in range(vmax+1)] for u in range(umax+1)] 
  s2y = [[0 for v in range(vmax+1)] for u in range(umax+1)] 
  s2z = [[0 for v in range(vmax+1)] for u in range(umax+1)] 
  #surface3
  s3x = [[0 for v in range(vmax+1)] for u in range(umax+1)] 
  s3y = [[0 for v in range(vmax+1)] for u in range(umax+1)] 
  s3z = [[0 for v in range(vmax+1)] for u in range(umax+1)] 
  #surface4
  s4x = [[0 for v in range(vmax+1)] for u in range(umax+1)] 
  s4y = [[0 for v in range(vmax+1)] for u in range(umax+1)] 
  s4z = [[0 for v in range(vmax+1)] for u in range(umax+1)] 

  splineSource1 = FindSource('SplineSource1')
  splineSource2 = FindSource('SplineSource2')
  splineSource3 = FindSource('SplineSource3')
  splineSource4 = FindSource('SplineSource4')

  pts1 = splineSource1.ParametricFunction.GetProperty("Points")
  pts2 = splineSource2.ParametricFunction.GetProperty("Points")
  pts3 = splineSource3.ParametricFunction.GetProperty("Points")
  pts4 = splineSource4.ParametricFunction.GetProperty("Points")

#==============================================
# linear interpolation based on 4 corner points 
#==============================================

  p1 = [pts1[0],pts1[1],pts1[2]] #start spline1
  p2 = [pts1[6],pts1[7],pts1[8]] #end   spline1
  p3 = [pts2[0],pts2[1],pts2[2]] #start spline2
  p4 = [pts2[6],pts2[7],pts2[8]] #end   spline2

  for uint in range(0, umax+1):
    px1 = p1[0] + (p2[0]-p1[0])*float(uint)/float(umax)
    px2 = p3[0] + (p4[0]-p3[0])*float(uint)/float(umax)
    py1 = p1[1] + (p2[1]-p1[1])*float(uint)/float(umax)
    py2 = p3[1] + (p4[1]-p3[1])*float(uint)/float(umax)
    pz1 = p1[2] + (p2[2]-p1[2])*float(uint)/float(umax)
    pz2 = p3[2] + (p4[2]-p3[2])*float(uint)/float(umax)
    dx = (px2-px1)/float(vmax)
    dy = (py2-py1)/float(vmax)
    dz = (pz2-pz1)/float(vmax)

    for vint in range(0, vmax+1):
      px = px1 + float(vint)*dx
      py = py1 + float(vint)*dy
      pz = pz1 + float(vint)*dz
      s1x[uint][vint] = px
      s1y[uint][vint] = py
      s1z[uint][vint] = pz
      # print uint,vint,px,py,pz
      # ps = PointSource(Center=[px,py,pz]) ; Show()

#========================================
# bezier interpolation based on 4 splines
#========================================

  #get points from spline point array
  p1 = [pts1[0],pts1[1],pts1[2]] #start spline1
  p2 = [pts1[3],pts1[4],pts1[5]] #mid
  p3 = [pts1[6],pts1[7],pts1[8]] #end

  p4 = [pts2[0],pts2[1],pts2[2]] #start spline2
  p5 = [pts2[3],pts2[4],pts2[5]] #mid
  p6 = [pts2[6],pts2[7],pts2[8]] #end

  p7 = [pts3[0],pts3[1],pts3[2]] #start spline3
  p8 = [pts3[3],pts3[4],pts3[5]] #mid
  p9 = [pts3[6],pts3[7],pts3[8]] #end

  p10 = [pts4[0],pts4[1],pts4[2]] #start spline4
  p11 = [pts4[3],pts4[4],pts4[5]] #mid
  p12 = [pts4[6],pts4[7],pts4[8]] #end

  #calc midpoint from spline points
  px1 = p1[0]*0.25 + p2[0]*0.5 + p3[0]*0.25
  py1 = p1[1]*0.25 + p2[1]*0.5 + p3[1]*0.25
  pz1 = p1[2]*0.25 + p2[2]*0.5 + p3[2]*0.25

  px2 = p4[0]*0.25 + p5[0]*0.5 + p6[0]*0.25
  py2 = p4[1]*0.25 + p5[1]*0.5 + p6[1]*0.25
  pz2 = p4[2]*0.25 + p5[2]*0.5 + p6[2]*0.25

  px3 = p7[0]*0.25 + p8[0]*0.5 + p9[0]*0.25
  py3 = p7[1]*0.25 + p8[1]*0.5 + p9[1]*0.25
  pz3 = p7[2]*0.25 + p8[2]*0.5 + p9[2]*0.25

  px4 = p10[0]*0.25 + p11[0]*0.5 + p12[0]*0.25
  py4 = p10[1]*0.25 + p11[1]*0.5 + p12[1]*0.25
  pz4 = p10[2]*0.25 + p11[2]*0.5 + p12[2]*0.25

  #move control point further out
  p2  = [ 3*p2[0] -2*px1, 3*p2[1] -2*py1, 3*p2[2] -2*pz1 ]
  p5  = [ 3*p5[0] -2*px2, 3*p5[1] -2*py2, 3*p5[2] -2*pz2 ]
  p8  = [ 3*p8[0] -2*px3, 3*p8[1] -2*py3, 3*p8[2] -2*pz3 ]
  p11 = [ 3*p11[0]-2*px4, 3*p11[1]-2*py4, 3*p11[2]-2*pz4 ]

  #========================
  # spline1 spline2 surface
  #========================

  for uint in range(0, umax+1):
    u = float(uint)/float(umax)
    px1 = p1[0]*(1-u)*(1-u) + 2.0*p2[0]*(1-u)*u + p3[0]*u*u
    py1 = p1[1]*(1-u)*(1-u) + 2.0*p2[1]*(1-u)*u + p3[1]*u*u
    pz1 = p1[2]*(1-u)*(1-u) + 2.0*p2[2]*(1-u)*u + p3[2]*u*u
    px2 = p4[0]*(1-u)*(1-u) + 2.0*p5[0]*(1-u)*u + p6[0]*u*u
    py2 = p4[1]*(1-u)*(1-u) + 2.0*p5[1]*(1-u)*u + p6[1]*u*u
    pz2 = p4[2]*(1-u)*(1-u) + 2.0*p5[2]*(1-u)*u + p6[2]*u*u
    dx = (px2-px1)/float(vmax)
    dy = (py2-py1)/float(vmax)
    dz = (pz2-pz1)/float(vmax)

    for vint in range(0, vmax+1):
      v  = float(vint)
      px = px1 + v*dx
      py = py1 + v*dy
      pz = pz1 + v*dz
      s2x[uint][vint] = px
      s2y[uint][vint] = py
      s2z[uint][vint] = pz
      # print uint,vint,px,py,pz
      # ps = PointSource(Center=[px,py,pz]);  Show()

  #========================
  # spline3 spline4 surface
  #========================

  for vint in range(0, vmax+1):
    v = float(vint)/float(vmax)
    px3 = p7[0]*(1-v)*(1-v)  + 2.0*p8[0]*(1-v)*v  + p9[0]*v*v
    py3 = p7[1]*(1-v)*(1-v)  + 2.0*p8[1]*(1-v)*v  + p9[1]*v*v
    pz3 = p7[2]*(1-v)*(1-v)  + 2.0*p8[2]*(1-v)*v  + p9[2]*v*v
    px4 = p10[0]*(1-v)*(1-v) + 2.0*p11[0]*(1-v)*v + p12[0]*v*v
    py4 = p10[1]*(1-v)*(1-v) + 2.0*p11[1]*(1-v)*v + p12[1]*v*v
    pz4 = p10[2]*(1-v)*(1-v) + 2.0*p11[2]*(1-v)*v + p12[2]*v*v
    dx = (px4 - px3)/float(umax)
    dy = (py4 - py3)/float(umax)
    dz = (pz4 - pz3)/float(umax)

    for uint in range(0, umax+1):
      u  = float(uint)
      px = px3 + u*dx
      py = py3 + u*dy
      pz = pz3 + u*dz
      s3x[uint][vint] = px
      s3y[uint][vint] = py
      s3z[uint][vint] = pz
      # print vint,uint,px,py,pz
      # ps = PointSource(Center=[px,py,pz]);  Show()

  #========================
  # surface4 = s2+s3-s1
  #========================

  for vint in range(0, vmax+1):
    for uint in range(0, umax+1):
      s4x[uint][vint] = s2x[uint][vint] + s3x[uint][vint] - s1x[uint][vint]
      s4y[uint][vint] = s2y[uint][vint] + s3y[uint][vint] - s1y[uint][vint]
      s4z[uint][vint] = s2z[uint][vint] + s3z[uint][vint] - s1z[uint][vint]

  #==========================
  # write surface to VTK file
  #==========================

  f = open("C:\zscratch\BESPOKE_CODES\SURFACE\surface.vtk","w")
  f.write("# vtk DataFile Version 3.0\n")
  f.write("vtk output\n")
  f.write("ASCII\n")
  f.write("DATASET STRUCTURED_GRID\n")
  f.write("DIMENSIONS %d %d %d \n" % (vmax+1,umax+1,1) )
  f.write("POINTS %d float\n" % npts)
  for uint in range(0, umax+1):
    for vint in range(0, vmax+1):
#     f.write("%f %f %f \n" % (s1x[uint][vint],s1y[uint][vint],s1z[uint][vint])) # corner points interp
#     f.write("%f %f %f \n" % (s2x[uint][vint],s2y[uint][vint],s2z[uint][vint])) # spline12 interpolate
#     f.write("%f %f %f \n" % (s3x[uint][vint],s3y[uint][vint],s3z[uint][vint])) # spline34 interpolate
      f.write("%f %f %f \n" % (s4x[uint][vint],s4y[uint][vint],s4z[uint][vint])) # bezier surface
  f.close()  

if True:
  polygon = LegacyVTKReader( FileNames=['C:\\zscratch\\BESPOKE_CODES\\SURFACE\\surface.vtk'] )
  DataRepresentation = Show()
  DataRepresentation.Representation = 'Surface With Edges'
  Render

