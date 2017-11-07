/* 
Copyright 2008 Guillaume ANCIAUX (guillaume.anciaux@epfl.ch)

This file is part of ParaViewHelper.

ParaViewHelper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ParaViewHelper is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ParaViewHelper.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static double * coordinates;
static int * connectivity;

inline int compareCoords(double * c1,double * c2){
  double tol=0.000001;
  double x1,x2,y1,y2,z1,z2; 

  x1 = c1[0];x2 = c2[0];
  if (fabs(x1-x2) > tol && x1 < x2) return -1;
  if (fabs(x1-x2) > tol && x2 < x1) return 1;
  
  y1 = c1[1];y2 = c2[1];
  if (fabs (y1-y2) > tol && y1 < y2) return -1;
  if (fabs (y1-y2) > tol && y2 < y1) return 1;

  z1 = c1[2];z2 = c2[2];
  if (fabs (z1-z2) > tol && z1 < z2) return -1;
  if (fabs (z1-z2) > tol && z2 < z1) return 1;
  
  return 0;
}


inline void AddNode(double * coords,int * nb_node, double x,double y,double z){
  coords[3*(*nb_node)+0]=x;
  coords[3*(*nb_node)+1]=y; 
  coords[3*(*nb_node)+2]=z;  
  (*nb_node)++;
}

inline void AddMedianNode(double *coords,int * nb_node,int n1,int n2){
  
  coords[3*(*nb_node)+0]= (coords[3*n1+0]+coords[3*n2+0])/2.;
  coords[3*(*nb_node)+1]= (coords[3*n1+1]+coords[3*n2+1])/2.;
  coords[3*(*nb_node)+2]= (coords[3*n1+2]+coords[3*n2+2])/2.;
  (*nb_node)++;
}

inline void AddElement(int * countelements,int * indexes,
		int n0,int n1,int n2,
		int n3,int n4,int n5,int n6,int n7,int n8,int n9){
  connectivity[10*(*countelements)+0]=indexes[n0]+1;
  connectivity[10*(*countelements)+1]=indexes[n1]+1;
  connectivity[10*(*countelements)+2]=indexes[n2]+1;
  connectivity[10*(*countelements)+3]=indexes[n3]+1;
  connectivity[10*(*countelements)+4]=indexes[n4]+1;
  connectivity[10*(*countelements)+5]=indexes[n5]+1;
  connectivity[10*(*countelements)+6]=indexes[n6]+1;
  connectivity[10*(*countelements)+7]=indexes[n7]+1;
  connectivity[10*(*countelements)+8]=indexes[n8]+1;
  connectivity[10*(*countelements)+9]=indexes[n9]+1;
  (*countelements)++;
}

inline int AddInSorted(double * node_coords,double * coords_array,int * index_array,int N){
  int low = 0;
  int high = N - 1;
  int mid=0;
  int i = 0;
  
  if (N == 0){
    coords_array[0] = node_coords[0];
    coords_array[1] = node_coords[1];
    coords_array[2] = node_coords[2];
    index_array[0] = -1;
    return 0;
  }

  while (low <= high)
    {
      mid = (low + high) / 2;
      switch (compareCoords(node_coords,coords_array+3*mid)){
      case -1: high = mid -1; break;
      case 1: low = mid + 1; break;
      case 0: return mid;break;//such a node as already been inserted
      default: fprintf(stderr,"AAAAAAAAARRRRRRRRRRG !! should not append\n");
      }
    }

  mid = (low + high) / 2;
  if (compareCoords(coords_array+3*mid,node_coords)==-1)
    ++mid;

  //shift all elements to insert at the right place the element we want
  
  double * c1 = coords_array+3*mid;
  double * c2;
  int * i1 = index_array+mid;
  int * i2;
  for (i = N-1; i >= mid ; --i){
    c1 = coords_array+3*i;
    c2 = c1+3;
    i1 = index_array+i;
    i2 = i1+1;
    c2[0] = c1[0];
    c2[1] = c1[1];
    c2[2] = c1[2];
    *i2 = * i1;
  }
  c1[0] = node_coords[0];
  c1[1] = node_coords[1];
  c1[2] = node_coords[2];
  *i1 = -1;

  return mid;
}


inline void AddNodesInGlobal(double * coords,int size,int * global_size,double * sorted_coordinates,int * sorted_indexes,int * indexes){
  //  int size2 = *global_size;
  int i;

  double c[3];
  
  for (i = 0 ; i < size ; ++i){
    //speedup search in sorted list to see if node already exists
    c[0] = coords[3*i];c[1] = coords[3*i+1];c[2] = coords[3*i+2];
    //    if (c[0] == .5 && c[1] == 1 && c[2] == 1) debug_flag = 1;

    int ind = AddInSorted(c,sorted_coordinates,sorted_indexes,*global_size);
    
    if (sorted_indexes[ind] == -1) {//node does not already exist

      coordinates[3*(*global_size)] = coords[3*i];
      coordinates[3*(*global_size)+1] = coords[3*i+1];
      coordinates[3*(*global_size)+2] = coords[3*i+2];
      indexes[i] = (*global_size);
      sorted_indexes[ind] = (*global_size);
      (*global_size)++;
    }
    else {//node already exist
      indexes[i] = sorted_indexes[ind];
    }
/*     fprintf(stderr,"rnodes: start\n");  */
/*     for (j = 0 ; j < *global_size ; ++j){ */
/*       fprintf(stderr,"(%e,%e,%e,%d)[%d]\n", */
/* 	      sorted_coordinates[3*j],sorted_coordinates[3*j+1],sorted_coordinates[3*j+2], */
/* 	      sorted_indexes[j],j);  */
/*     } */
/*     fprintf(stderr,"end\n");  */
  }
}

void write_mesh_UNV(const char* outputfile,int nodes,int elements){
    int i;
    FILE *op = fopen(outputfile, "w");

    //nodal values (not used at present time)
    int exp_coord_sys_num = 0;
    int disp_coord_sys_num = 0;
    int color = 0;

    //element values 
    int fe_descriptor_id = 111;
    int phys_prop_tab_num = 2;
    int mat_prop_tab_num = 1;
    int n_nodes = 4;
    
    //start of node block
    fprintf(op,"    -1\n  2411\n");
    
    for (i=0;i<nodes;i++){
      fprintf(op,"         %d         %d         %d         %d\n",
	      i,exp_coord_sys_num,disp_coord_sys_num,color);
      fprintf(op,"  %.15e  %.15e  %.15e\n",
                    coordinates[3*i+0],
                    coordinates[3*i+1],
                    coordinates[3*i+2]);
    }
    //end of node block
    fprintf(op,"    -1\n");
    //start of element block
    fprintf(op,"    -1\n  2412\n");
    color = 7;
    for (i=0;i<elements;i++){
      fprintf(op,"         %d         %d         %d         %d         %d         %d\n",
	      i,fe_descriptor_id,phys_prop_tab_num,mat_prop_tab_num,color,n_nodes);
      fprintf(op,"         %d         %d         %d         %d\n",
          connectivity[10*i+0]-1,connectivity[10*i+1]-1,connectivity[10*i+2]-1,connectivity[10*i+3]-1);
    }
    //end of element block
    fprintf(op,"-1\n");

    fclose(op);
}


void BuildGridMesh(double ** pos,int ** conn, 
		   int * nodes,int * elements,
		   double L, double D, double W,
		   double NL,double ND, double NW)
{
  int i,j,k;
  int countnodes,countelements;
  int nstart,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9;
  double lcube,dcube,wcube;
  double x0,y0,z0;
  //  double tol = 0.00001;
  
  lcube = L/NL;
  dcube = D/ND;
  wcube = W/NW;
  *nodes = 35*NL*ND*NW; /* 35 is the number of nodes per unit cube; node is overestimated but identical nodes will be removed later on */
  *elements = 12*NL*ND*NW;
  
  *pos = malloc(*nodes*3*sizeof(double));
  coordinates = *pos;
  double * sorted_coordinates = malloc(*nodes*3*sizeof(double));
  int * sorted_indexes = malloc(*nodes*sizeof(int));
  *conn = malloc(*elements*10*sizeof(int));
  connectivity = *conn;


  //  int global_size=0;

  countnodes=0;
  countelements=0;

  for (i=0;i<NW;i++)
    for (j=0;j<ND;j++)
      for (k=0;k<NL;k++)
	{

	  if (countelements%1000 == 0)
	    fprintf(stderr,"building element %d/%d\n",countelements,*elements);
	  x0 = lcube*k;
	  y0 = dcube*j;
	  z0 = wcube*i;

	  double temp_coords[3*35];
	  int indexes[35];
	  int cpt=0;
	  /* Base corner nodes */
	  nstart = countnodes+1;
	  AddNode(temp_coords,&cpt,x0,y0,z0);
	  AddNode(temp_coords,&cpt,x0+lcube,y0,z0);
	  AddNode(temp_coords,&cpt,x0+lcube,y0+dcube,z0);
	  AddNode(temp_coords,&cpt,x0,y0+dcube,z0);
	  /* Top corner nodes */
	  AddNode(temp_coords,&cpt,x0,y0,z0+wcube);
	  AddNode(temp_coords,&cpt,x0+lcube,y0,z0+wcube);
	  AddNode(temp_coords,&cpt,x0+lcube,y0+dcube,z0+wcube);
	  AddNode(temp_coords,&cpt,x0,y0+dcube,z0+wcube);
	  
	  /* Mid node */
	  n1=0;n2=6;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
 	  /* Mid nodes on the faces */
	  n1=0;n2=1;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=1;n2=2;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=2;n2=3;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=0;n2=3;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=4;n2=5;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=5;n2=6;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=6;n2=7;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=7;n2=4;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=0;n2=4;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=1;n2=5;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=2;n2=6;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=3;n2=7;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  
	  /* Mid nodes on the diagonals of the faces */
	  n1=4;n2=1;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=6;n2=1;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=7;n2=2;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=7;n2=0;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=3;n2=1;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=7;n2=5;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  
	  /* Joining the corner nodes to the central mid node n9 */
	  n1=0;n2=8;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=1;n2=8;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=2;n2=8;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=3;n2=8;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=4;n2=8;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=5;n2=8;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=6;n2=8;
	  AddMedianNode(temp_coords,&cpt,n1,n2);
	  n1=7;n2=8;
	  AddMedianNode(temp_coords,&cpt,n1,n2);

	  //memcpy(coordinates+(nstart-1)*3,temp_coords,3*nb_nodes*sizeof(double));
	  AddNodesInGlobal(temp_coords,cpt,&countnodes,sorted_coordinates,sorted_indexes,indexes);
	  //	  countnodes += cpt;
	  
	  /* CONNECTIVITY */
	  n0=0;n1=4;n2=1;n3=8;n4=17;n5=21;n6=9;n7=27;n8=31;n9=28;
	  AddElement(&countelements,indexes,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9);
	  
	  n0=1;n1=4;n2=5;n3=8;n4=21;n5=13;n6=18;n7=28;n8=31;n9=32;
	  AddElement(&countelements,indexes,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9);

	  n0=1;n1=5;n2=6;n3=8;n4=18;n5=14;n6=22;n7=28;n8=32;n9=33;
	  AddElement(&countelements,indexes,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9);

	  n0=1;n1=6;n2=2;n3=8;n4=22;n5=19;n6=10;n7=28;n8=33;n9=29;
	  AddElement(&countelements,indexes,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9);

	  n0=2;n1=6;n2=7;n3=8;n4=19;n5=15;n6=23;n7=29;n8=33;n9=34;
	  AddElement(&countelements,indexes,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9);

	  n0=2;n1=7;n2=3;n3=8;n4=23;n5=20;n6=11;n7=29;n8=34;n9=30;
	  AddElement(&countelements,indexes,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9);

	  n0=3;n1=7;n2=0;n3=8;n4=20;n5=24;n6=12;n7=30;n8=34;n9=27;
	  AddElement(&countelements,indexes,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9);

	  n0=7;n1=4;n2=0;n3=8;n4=16;n5=17;n6=24;n7=34;n8=31;n9=27;
	  AddElement(&countelements,indexes,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9);

	  n0=4;n1=7;n2=5;n3=8;n4=16;n5=26;n6=13;n7=31;n8=34;n9=32;
	  AddElement(&countelements,indexes,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9);

	  n0=7;n1=6;n2=5;n3=8;n4=15;n5=14;n6=26;n7=34;n8=33;n9=32;
	  AddElement(&countelements,indexes,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9);

	  n0=0;n1=1;n2=3;n3=8;n4=9;n5=25;n6=12;n7=27;n8=28;n9=30;
	  AddElement(&countelements,indexes,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9);

	  n0=1;n1=2;n2=3;n3=8;n4=10;n5=11;n6=25;n7=28;n8=29;n9=30;
	  AddElement(&countelements,indexes,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9);
	}
  
  *nodes = countnodes;
  //  write_mesh_UNV("cube.unv",*nodes,*elements);
}
