#ifndef WRITE_MATLAB_ARRAY_H
#define WRITE_MATLAB_ARRAY_H

// ------------------------------------------------------------------------------------
// Save an array to a matlab file.
// ------------------------------------------------------------------------------------

int writeMatlabArray(FILE *matlabFile, Real *u_p, const char *name, int nd1a, int nd1b, int nd2a, int nd2b)
{
	#define u(i1, i2) u_p[(i1 - nd1a) + nd1 * (i2 - nd2a)]
	const int nd1 = nd1b-nd1a+1;
	const int nd2 = nd2b-nd2a+1;

	const int numPerLine=8; // number of entries per line

	fprintf(matlabFile, "%s=zeros(%d,%d);\n", name, nd1, nd2);
  int count = 0;
  for(int i2 = nd2a; i2 <= nd2b; i2++)
		{
			for( int i1 = nd1a; i1 <= nd1b; i1++ )
			{
				fprintf(matlabFile, "%s(%3d,%3d)=%12.5e; ", name, i1-nd1a+1,i2-nd2a+1,u(i1,i2));
				if( count % numPerLine == numPerLine-1 )
					fprintf(matlabFile,"\n"); // new line
				count++;
			}
		}
	fprintf(matlabFile,"\n");
return 0;
}

#endif