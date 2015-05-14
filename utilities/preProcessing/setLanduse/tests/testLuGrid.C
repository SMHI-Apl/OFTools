#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include "Matrix.H"
#include "luGrid.H"
main()
{
  int x, y;
  ifstream fid("landuse_grid.txt");
  ofstream fut("landuse_ascii_ut.txt");
  luGrid m(fid);
 
  cout << "Grid size is: " << endl;
  cout << "Antal rader: " << m.numRow()<< " Antal kolumner: "<< m.numCol()<< endl;
  cout << "cellSize is: "<<m.getCellSize()<<endl;
  
  for(int i=1;i<21;i++){
    for(int j=1;j<21;j++)	
      cout<< m(i,j) <<" ";
    cout<<endl;
  }

  cout<< "Code in -499, -499 is: "<<m.getCode(-499,-499)<<endl;
 
  m.write(fut);

  return 0;
}
