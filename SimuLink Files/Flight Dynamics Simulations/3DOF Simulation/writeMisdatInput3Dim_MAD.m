function [status] = writeMisdatInput3Dim_MAD(Sdim,alt)
  %Write MISDAT for005.dat (assumes 3 dimensions)
  
  fileName = 'for005.dat';
  
  fid = fopen(fileName,'w');
  
  fprintf(fid,'CASEID TEST CASE FOR KRIGING\n');
  fprintf(fid,'$FLTCON\n');
  fprintf(fid,' NALPHA = 1.,\n');
  fprintf(fid,' ALPHA = %6.4f,\n',Sdim(1,2));
  fprintf(fid,' NMACH = 1.,\n');
  fprintf(fid,' MACH = %6.4f,\n',Sdim(1,1));
  fprintf(fid,' ALT = %6.1f,\n',alt);
  fprintf(fid,'$END\n');
  fprintf(fid,'$REFQ\n');
  fprintf(fid,' XCG = 2.00,\n');
  fprintf(fid,'$END\n');
  fprintf(fid,'$AXIBOD\n');
  fprintf(fid,' TNOSE = 1.,\n');
  fprintf(fid,' LNOSE = 0.9,\n');
  fprintf(fid,' DNOSE = 0.3,\n');
  fprintf(fid,' LCENTR = 3.03,\n');
  fprintf(fid,' DCENTR = 0.30,\n');
  fprintf(fid,' TAFT = 0.,\n');
  fprintf(fid,' LAFT = 0.12,\n');
  fprintf(fid,' DAFT = 0.22,\n');
  fprintf(fid,' DEXIT = 0.,\n');
  fprintf(fid,'$END\n');
  fprintf(fid,'$FINSET1\n');
  fprintf(fid,' XLE = 3.48,\n');
  fprintf(fid,' NPANEL = 4.,\n');
  fprintf(fid,' PHIF = 45., 135., 225., 315.,\n');
  fprintf(fid,' SWEEP = 50.44,\n');
  fprintf(fid,' STA = 0.,\n');
  fprintf(fid,' SSPAN = 0.0, 0.19,\n');
  fprintf(fid,' CHORD = 0.45, 0.225,\n');
  fprintf(fid,' ZUPPER = 0.025, 0.025,\n');
  fprintf(fid,' ZLOWER = 0.025, 0.025,\n');
  fprintf(fid,'$END\n');
  fprintf(fid,'$DEFLCT\n');
  fprintf(fid,' DELTA1 = %6.4f, %6.4f, %6.4f, %6.4f,\n',-Sdim(1,3),-Sdim(1,3),Sdim(1,3),Sdim(1,3));
  fprintf(fid,'$END\n');
  fprintf(fid,'SAVE\n');
  fprintf(fid,'DIM M\n');
  fprintf(fid,'NEXT CASE\n');
  
  for iloop = 2:length(Sdim)
    fprintf(fid,'$FLTCON\n');
    fprintf(fid,' NALPHA = 1.,\n');
    fprintf(fid,' ALPHA = %6.4f,\n',Sdim(iloop,2));
    fprintf(fid,' NMACH = 1.,\n');
    fprintf(fid,' MACH = %6.4f,\n',Sdim(iloop,1));
    fprintf(fid,'$END\n');
    fprintf(fid,'$DEFLCT\n');
    fprintf(fid,' DELTA1 = %6.4f, %6.4f, %6.4f, %6.4f,\n',-Sdim(iloop,3),-Sdim(iloop,3),Sdim(iloop,3),Sdim(iloop,3));
    fprintf(fid,'$END\n');
    fprintf(fid,'SAVE\n');
    fprintf(fid,'NEXT CASE\n');
  end
  
  fclose(fid);
  
  status = 1;