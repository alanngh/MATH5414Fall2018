/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "parAlmond.hpp"



namespace parAlmond {

// a < b
bool CompareNBS(dlong Na, dlong oNa, dlong Ia, dlong Nb, dlong oNb, dlong Ib){
	if (Na < Nb) 	return true;
	if (Na > Nb) 	return false;
	if (oNa > oNb)	return true;
	if (oNa < oNb)	return false;
	if (Ia > Ib)	return true;
	if (Ia < Ib)	return false;
	
	return false; 
}



void formAggregates(parCSR *A, parCSR *C,
                     hlong* FineToCoarse,
                     hlong* globalAggStarts){

  int rank, size;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &size);

  const dlong N   = C->Nrows;
  const dlong M   = C->Ncols;
  

  printf("----------------------\n");
  printf(" rank = %d   N =  %d M = %d \n",rank,N,M);
  printf("----------------------\n");

	int Rprint = 1;

  int   *rands = (int *)     calloc(M, sizeof(int));  // initialized to 0s
  int   *states = (int *)    calloc(M, sizeof(int));     // initialized to 0s
  int   *indices = (int *)   calloc(M, sizeof(int));     // initialized to 0s

  // hlong *globalRowStarts = A->globalRowStarts;

  for(dlong i=0; i<N; i++)
	rands[i] = C->diag->rowStarts[i+1] - C->diag->rowStarts[i] + C->offd->rowStarts[i+1] - C->offd->rowStarts[i] ;
 
  //sort indices
  for (int i=0;i<N;i++) 	
	indices[i] = i;
		  
  for (int i=0;i<N;i++){
	for (int j=i+1;j<N;j++){
		dlong oNa = C->offd->rowStarts[indices[i]+1]-C->offd->rowStarts[indices[i]];
		dlong oNb = C->offd->rowStarts[indices[j]+1]-C->offd->rowStarts[indices[j]];
		if (CompareNBS(rands[indices[i]],oNa,indices[i],rands[indices[j]],oNb,indices[j])){
			//printf("\n (%d , %d , %d)  < (%d , %d , %d)",rands[i],oNa,i,rands[j],oNb,j);
			int temp = indices[i];
			indices[i] = indices[j];
			indices[j] =  temp;			
		}		
	} 
  }
   
  bool done = true;
  hlong Cnt=0, lCnt, oldCnt;
  
  
  while(done){
	
	lCnt = std::count(states,states+N,0);
	MPI_Allreduce(&lCnt,&oldCnt,1,MPI_HLONG,MPI_SUM,A->comm);
  
	for (int i=0;i<N;i++){
		// verifica que la vecindad esta disponible
		if (states[indices[i]] == 0 ){
			bool ok = true;
			if (ok){ // no local
				for (dlong j = C->offd->rowStarts[indices[i]]; j<C->offd->rowStarts[indices[i]+1];j++){
					if (states[C->offd->cols[j]] != 0){
						if (states[C->offd->cols[j]] == 1)
							states[indices[i]] = -1;
						ok = false;
						//break;	revisa todos los vecinos para saber si esta conectado a un nodo raiz 			
					}				
				}			
			}
			if (ok){// local
				for (dlong j = C->diag->rowStarts[indices[i]]; j<C->diag->rowStarts[indices[i]+1];j++){
					if (states[C->diag->cols[j]] != 0){
						ok = false;
						break;				
					}				
				}
			}		
			if (ok){
				for (dlong j = C->diag->rowStarts[indices[i]]; j<C->diag->rowStarts[indices[i]+1];j++){
					if (C->diag->cols[j] ==  indices[i])		
						states[C->diag->cols[j]] = 1;				
					else 
						states[C->diag->cols[j]] = -1;				
				}
				
				//printf("\n----------\n RAnk = %d Aggregates created at node %d -\n-----------\n",rank,indices[i]);									
				break;				
			}		
		}						
	}
	
	for (int n = N ; n<M ; n++) states[n]=0;
	ogsGatherScatter(states,ogsInt,ogsAdd,A->ogs);		

	
	lCnt = std::count(states,states+N,0);
	MPI_Allreduce(&lCnt,&Cnt,1,MPI_HLONG,MPI_SUM,A->comm);
	
	if(oldCnt == Cnt)	done = false;
 }   
  
  // printf("\n rank %d  tiene %d/%d sin asignar\n",rank,lCnt,Cnt);
   
   dlong numAggs = 0;
   dlong *gNumAggs = (dlong *) calloc(size,sizeof(dlong));
   
   for ( dlong i=0; i<N ; i++)
		if (states[i]==1)	numAggs++;	
	  
   MPI_Allgather(&numAggs,1,MPI_DLONG,gNumAggs,1,MPI_DLONG,A->comm);
   
   globalAggStarts[0] = 0;
   for (int r=0;r<size;r++)
		globalAggStarts[r+1] = globalAggStarts[r] + gNumAggs[r];
				
	numAggs = 0;
	
	for (int i=0;i<N ; i++){
		if (states[indices[i]]==1)	
			FineToCoarse[indices[i]] = globalAggStarts[rank] + numAggs++;
		else
			FineToCoarse[indices[i]] = -1;
	}
   
   for (dlong i=N ; i<M ; i++)	FineToCoarse[i] = 0;
   ogsGatherScatter(FineToCoarse,ogsHlong,ogsAdd,A->ogs);

/*
	if (rank==1){
	printf("\n initial flags...\n");
	for (int i=0;i<N;i++)	
		printf("%d  ",FineToCoarse[i]);
	}
*/
     
   // cosntruct local aggregates
   for (int i=0;i<N;i++){
		if (states[i] == 1 ){
			for (dlong j = C->diag->rowStarts[i]; j<C->diag->rowStarts[i+1];j++){
				FineToCoarse[C->diag->cols[j]] = FineToCoarse[i];
			}
		}				
	}
	
/*	if (rank==1){
	printf("\n local aggs...\n");
	for (int i=0;i<N;i++)	
		printf("%d  ",FineToCoarse[i]);
	}
  */ 
   
   for (dlong i=N ; i<M ; i++)	FineToCoarse[i] = 0;
   ogsGatherScatter(FineToCoarse,ogsHlong,ogsAdd,A->ogs);
   
   // cosntruct non-local aggregates
   for (int i=0;i<N;i++){
		if (FineToCoarse[indices[i]] < 0 ){
			for (dlong j = C->offd->rowStarts[indices[i]] ; j<C->offd->rowStarts[indices[i]+1] ; j++){
				if (FineToCoarse[C->offd->cols[j]] > 0){ // empaeja al primer vecino 
					FineToCoarse[indices[i]] = FineToCoarse[C->offd->cols[j]];
					break;
				}
			}
		}				
	}
   
   
   
/*	if (rank==1){
	printf("\n non-local aggs..\n");
	for (int i=0;i<N;i++)	
		printf("%d  ",FineToCoarse[i]);
	}*/
   
   
   for (dlong i=N ; i<M ; i++)	FineToCoarse[i] = 0;
   ogsGatherScatter(FineToCoarse,ogsHlong,ogsAdd,A->ogs);
   
   int mAgg = 0;
   
   for (int i=0 ; i<N ; i++)
		if (FineToCoarse[i]<0)		mAgg++;
		
//	printf("\n rank =%d   with %d missing nodes !!!  \n",rank,mAgg);
	   
   if (mAgg>0){
	  bool ok = false;
	  for (int i=0 ; i<N ; i++){
		  if (FineToCoarse[i]< 0 ){
			for (dlong j = C->diag->rowStarts[i]; j<C->diag->rowStarts[i+1];j++){
				if (FineToCoarse[C->diag->cols[j]] > -1){
						FineToCoarse[i] = FineToCoarse[C->diag->cols[j]];
						ok = true;
						break;
				}							
			}
			if (!ok){
				for (dlong j = C->offd->rowStarts[i]; j<C->offd->rowStarts[i+1];j++){
					if (FineToCoarse[C->diag->cols[j]] > -1){
						FineToCoarse[i] = FineToCoarse[C->diag->cols[j]];
						ok = true;
						break;
					}							
				}					
			}
			if (!ok){
				printf("\n Rank %d missing node %d .... fck of !!!",rank,i);	    
			}
		}
	  }
   }
   
   /*
   if (rank==1){
	printf("\n final aggs..\n");
	for (int i=0;i<N;i++)	
		printf("%d  ",FineToCoarse[i]);
	}
   */
      
   for (dlong i=N ; i<M ; i++)	FineToCoarse[i] = 0;   
   ogsGatherScatter(FineToCoarse,ogsHlong,ogsAdd,A->ogs);
	
	
	
	free(states);
	delete C;
	
}


void formAggregates2(parCSR *A, parCSR *C,
                     hlong* FineToCoarse,
                     hlong* globalAggStarts){

  int rank, size;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &size);

  const dlong N   = C->Nrows;
  const dlong M   = C->Ncols;
  const dlong diagNNZ = C->diag->nnz;
  const dlong offdNNZ = C->offd->nnz;

	printf("----------------------\n");
	printf(" rank = %d   N =  %d M = %d \n",rank,N,M);
	printf("----------------------\n");
	int Rprint = -1;



  dfloat *rands = (dfloat *) calloc(M, sizeof(dfloat));
  int   *states = (int *)    calloc(M, sizeof(int));

  dfloat *Tr = (dfloat *) calloc(M, sizeof(dfloat));
  int    *Ts = (int *)    calloc(M, sizeof(int));
  hlong  *Ti = (hlong *)  calloc(M, sizeof(hlong));
  hlong  *Tc = (hlong *)  calloc(M, sizeof(hlong));

  hlong *globalRowStarts = A->globalRowStarts;

  for(dlong i=0; i<N; i++)
	rands[i] = (dfloat) drand48();

  // add the number of non-zeros in each column
  int *colCnt = (int *) calloc(M,sizeof(int));
  for(dlong i=0; i<diagNNZ; i++)
    colCnt[C->diag->cols[i]]++;

  for(dlong i=0; i<offdNNZ; i++)
    colCnt[C->offd->cols[i]]++;

  //gs for total column counts
  ogsGatherScatter(colCnt, ogsInt, ogsAdd, A->ogs);

  //add random pertubation
  for(int i=0;i<N;++i)
    rands[i] += colCnt[i];

  //gs to fill halo region
  ogsGatherScatter(rands, ogsDfloat, ogsAdd, A->ogs);

  hlong done = 0;
  while(!done){
    // first neighbours
    // #pragma omp parallel for
    for(dlong i=0; i<N; i++){

      int smax = states[i];
      dfloat rmax = rands[i];
      hlong imax = i + globalRowStarts[rank];

      if(smax != 1){
        //local entries
        for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
          const dlong col = C->diag->cols[jj];
          if (col==i) continue;
          if(customLess(smax, rmax, imax, states[col], rands[col], col + globalRowStarts[rank])){
            smax = states[col];
            rmax = rands[col];
            imax = col + globalRowStarts[rank];
          }
        }
        //nonlocal entries
        for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
          const dlong col = C->offd->cols[jj];
          if(customLess(smax, rmax, imax, states[col], rands[col], A->colMap[col])) {
            smax = states[col];
            rmax = rands[col];
            imax = A->colMap[col];
          }
        }
      }
      Ts[i] = smax;
      Tr[i] = rmax;
      Ti[i] = imax;
    }

    //share results
    for (dlong n=N;n<M;n++) {
      Tr[n] = 0.;
      Ts[n] = 0;
      Ti[n] = 0;
    }
    ogsGatherScatter(Tr, ogsDfloat, ogsAdd, A->ogs);
    ogsGatherScatter(Ts, ogsInt,    ogsAdd, A->ogs);
    ogsGatherScatter(Ti, ogsHlong,  ogsAdd, A->ogs);

    // second neighbours
    // #pragma omp parallel for
    for(dlong i=0; i<N; i++){
      int    smax = Ts[i];
      dfloat rmax = Tr[i];
      hlong  imax = Ti[i];

      //local entries
      for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
        const dlong col = C->diag->cols[jj];
        if (col==i) continue;
        if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
          smax = Ts[col];
          rmax = Tr[col];
          imax = Ti[col];
        }
      }
      //nonlocal entries
      for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
        const dlong col = C->offd->cols[jj];
        if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
          smax = Ts[col];
          rmax = Tr[col];
          imax = Ti[col];
        }
      }

      // if I am the strongest among all the 1 and 2 ring neighbours
      // I am an MIS node
      if((states[i] == 0) && (imax == (i + globalRowStarts[rank])))
        states[i] = 1;

      // if there is an MIS node within distance 2, I am removed
      if((states[i] == 0) && (smax == 1))
        states[i] = -1;
    }

    //share results
    for (dlong n=N;n<M;n++) states[n] = 0;
    ogsGatherScatter(states, ogsInt, ogsAdd, A->ogs);

    // if number of undecided nodes = 0, algorithm terminates
    hlong cnt = std::count(states, states+N, 0);
    MPI_Allreduce(&cnt,&done,1,MPI_HLONG, MPI_SUM,A->comm);
    done = (done == 0) ? 1 : 0;
  }
  
  /*
  if (rank ==Rprint){
	  printf("\n   root nodes \n");
	for (dlong i =0 ; i<M ; i++) 	printf("\n %d    %d  ",i,states[i]);
  }
  */
  
  dlong numAggs = 0;
  dlong *gNumAggs = (dlong *) calloc(size,sizeof(dlong));

  // count the coarse nodes/aggregates
  for(dlong i=0; i<N; i++)
    if(states[i] == 1) numAggs++;

  MPI_Allgather(&numAggs,1,MPI_DLONG,gNumAggs,1,MPI_DLONG,A->comm);

  globalAggStarts[0] = 0;
  for (int r=0;r<size;r++)
    globalAggStarts[r+1] = globalAggStarts[r] + gNumAggs[r];

  numAggs = 0;
  // enumerate the coarse nodes/aggregates
  for(dlong i=0; i<N; i++) {
    if(states[i] == 1) {
      FineToCoarse[i] = globalAggStarts[rank] + numAggs++;
    } else {
      FineToCoarse[i] = -1;
    }
  }
  for(dlong i=N; i<M; i++) FineToCoarse[i] = 0;

  //share the initial aggregate flags
  ogsGatherScatter(FineToCoarse, ogsHlong, ogsAdd, A->ogs);
  
  /*
   if (rank == Rprint){
	   printf("\n initial flags \n");
	for (dlong i =0 ; i<M ; i++) 	printf("\n %d    %d  ",i,FineToCoarse[i]);
  }*/
  
  // form the aggregates
  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    int   smax = states[i];
    dfloat rmax = rands[i];
    hlong  imax = i + globalRowStarts[rank];
    hlong  cmax = FineToCoarse[i];

    if(smax != 1){
      //local entries
      for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
        const dlong col = C->diag->cols[jj];
        if (col==i) continue;
        if(customLess(smax, rmax, imax, states[col], rands[col], col + globalRowStarts[rank])){
          smax = states[col];
          rmax = rands[col];
          imax = col + globalRowStarts[rank];
          cmax = FineToCoarse[col];
        }
      }
      //nonlocal entries
      for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
        const dlong col = C->offd->cols[jj];
        if(customLess(smax, rmax, imax, states[col], rands[col], A->colMap[col])){
          smax = states[col];
          rmax = rands[col];
          imax = A->colMap[col];
          cmax = FineToCoarse[col];
        }
      }
    }
    Ts[i] = smax;
    Tr[i] = rmax;
    Ti[i] = imax;
    Tc[i] = cmax;

    if((states[i] == -1) && (smax == 1) && (cmax > -1))
      FineToCoarse[i] = cmax;
  }

  //share results
  for (dlong n=N;n<M;n++) {
    FineToCoarse[n] = 0;
    Tr[n] = 0.;
    Ts[n] = 0;
    Ti[n] = 0;
    Tc[n] = 0;
  }
  
  /* if (rank == Rprint){
	   printf("\n first phase (before GS) \n");
	for (dlong i =0 ; i<M ; i++) 	printf("\n %d    %d  ",i,FineToCoarse[i]);
  }*/
   
  ogsGatherScatter(FineToCoarse, ogsHlong,  ogsAdd, A->ogs);
  ogsGatherScatter(Tr,     ogsDfloat, ogsAdd, A->ogs);
  ogsGatherScatter(Ts,     ogsInt,    ogsAdd, A->ogs);
  ogsGatherScatter(Ti,     ogsHlong,  ogsAdd, A->ogs);
  ogsGatherScatter(Tc,     ogsHlong,  ogsAdd, A->ogs);
  
  /*
  
  if (rank == Rprint){
	  printf("\n firs phase (after GS) \n");
	for (dlong i =0 ; i<M ; i++) 	printf("\n %d    %d  ",i,FineToCoarse[i]);
  }*/
  
  
  
  // second neighbours
  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    int    smax = Ts[i];
    dfloat rmax = Tr[i];
    hlong  imax = Ti[i];
    hlong  cmax = Tc[i];

    //local entries
    for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
      const dlong col = C->diag->cols[jj];
      if (col==i) continue;
      if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
        smax = Ts[col];
        rmax = Tr[col];
        imax = Ti[col];
        cmax = Tc[col];
      }
    }
    //nonlocal entries
    for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
      const dlong col = C->offd->cols[jj];
      if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
        smax = Ts[col];
        rmax = Tr[col];
        imax = Ti[col];
        cmax = Tc[col];
      }
    }

    if((states[i] == -1) && (smax == 1) && (cmax > -1))
      FineToCoarse[i] = cmax;
  }
 
 
  /*if (rank == Rprint){
	  printf("\n second phase (before GS) \n");
	for (dlong i =0 ; i<M ; i++) 	printf("\n %d    %d  ",i,FineToCoarse[i]);
  }*/
  
 
  //share results
  for (dlong n=N;n<M;n++) FineToCoarse[n] = 0;
  ogsGatherScatter(FineToCoarse, ogsHlong,  ogsAdd, A->ogs);

 /*if (rank == Rprint){
	  printf("\n second phase (after GS) \n");
	for (dlong i =0 ; i<M ; i++) 	printf("\n %d    %d  ",i,FineToCoarse[i]);
  }*/
  


  free(rands);
  free(states);
  free(Tr);
  free(Ts);
  free(Ti);
  free(Tc);

  delete C;
}


} //namespace parAlmond
