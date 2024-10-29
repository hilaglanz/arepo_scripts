fin  = open( "tracer.dat", "rb" )
fout = open( "tracer_serial.dat", "wb" )

fin.seek( 4, 0 )
NTracer, = struct.unpack( "i", fin.read(4) )

fin.seek( 8, 1 )
Masses = np.fromfile( fin, count=NTracer, dtype='f8' )

offset   = NTracer * 8 + 20
tstpsize = 16 + 6 * 4 * NTracer

print( "NTracer: ", NTracer )
print( "Total Mass: ", Masses.sum()/msol )

Times = np.zeros( 10**6, dtype='f4' )

# find number of timesteps
NTimesteps = 0
fin.seek( offset + NTimesteps * tstpsize + 4, 0 )
s = fin.read(8)

start = time.time()

while len(s) > 0:
  Time, = struct.unpack( "d", s )
  Times[NTimesteps] = Time
  NTimesteps += 1
  fin.seek( offset + NTimesteps * tstpsize + 4, 0 )
  s = fin.read(8)
  
  if NTimesteps % 10000 == 0:
    print( "Timesteps currently: %d" % NTimesteps )

print( "Total number of timesteps: ", NTimesteps )
print( "Time runs from ", Times[0], " to ", Times[NTimesteps-1] )

print( "Getting timesteps took %ds." % (time.time() - start) )

fout.write( struct.pack( "iiiii", NTracer, 1, NTracer, NTimesteps, 6 ) )
Masses.tofile( fout )
Times[:NTimesteps].tofile( fout )

# 1GB buffer
ReadChunk = 2 * 10**11 // (6 * NTimesteps * 8)

print( "Reading ", ReadChunk, " tracers per chunk." )

ReadCount = 0
while ReadCount < NTracer:
  ReadChunk = min( ReadChunk, NTracer - ReadCount )
  
  print( ReadChunk )
  
  data = np.zeros( (ReadChunk,NTimesteps,6), dtype='f4' )
  
  start = time.time()
  
  print( "Reading chunk size %d starting from %d." % (ReadChunk, ReadCount) )
  for itstp in range( NTimesteps ):
    if itstp % 10000 == 0:
      print( "Now reading timestep %d/%d." % (itstp,NTimesteps) )
    
    for ival in range(6):
      fin.seek( offset + itstp * tstpsize + 12 + NTracer * ival * 4 + ReadCount * 4, 0 )
      chunk = np.fromfile( fin, dtype='f4', count=ReadChunk )
      data[:,itstp,ival] = chunk
  
  print( "Writing chunk size %d starting from %d." % (ReadChunk, ReadCount) )
  for itracer in range( ReadChunk ):
    print( "tracer %6d: rho=%g - %g" % (itracer+ReadCount, data[itracer,:,3].min(), data[itracer,:,3].max()) )
    data[itracer,:,:].T.tofile( fout )

  print( "Writing done, reading and writing of chunk took %ds." % (time.time() - start) )
  
  ReadCount += ReadChunk
  
fin.close()
fout.close()
